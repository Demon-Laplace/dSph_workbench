from multiprocessing import get_context
import argparse
import ctypes
import os
import resource
import time
from basefunc import DataProcessor, GalaxySimulation
from snapshot_context import prepare_snapshot_context
from snapshot_metrics import compute_snapshot_summary
from variable import core_radius_dict, folder_path, gizmo_config, orbit_file
from pathlib import Path
import sys
import csv
import gc
import tempfile
from contextlib import contextmanager

try:
    import fcntl
except ImportError:
    fcntl = None

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

ELINFO_COLUMNS = [
    'star_mass', 'star_half_mass',
    'hotgas_mass', 'hotgas_half_mass',
    'coldgas_mass', 'coldgas_half_mass',
    'mw_mass_r', 'eps', 'pa', 'gas_density',
    'vr', 'vtheta', 'vphi', 'distance_gal', 'rhalf',
    'rhalf_circularized', 'rhalf_circular',
    'shape_center_x_kpc', 'shape_center_y_kpc',
    'distance', 'age', 'sigma', 'sigma_re_circular',
    'sigma_re_elliptical', 'sigma_fixed_500pc', 'sigma_re_nstar',
    'sigma_gradient_kms_per_kpc', 'tsigma',
    'tsigma_xyz',
    'sigma_x', 'sigma_y', 'sigma_z',
    'pmra', 'pmdec', 'cold_gas_center_ra', 'cold_gas_center_dec',
    'numsp',
]
NUMSP_COLUMN = 'numsp'
NUMSP_INDEX = ELINFO_COLUMNS.index(NUMSP_COLUMN)

try:
    LIBC = ctypes.CDLL("libc.so.6")
except OSError:
    LIBC = None

try:
    from threadpoolctl import threadpool_limits
except ImportError:
    threadpool_limits = None

THREADPOOL_CONTROLLER = None


def build_consecutive_blocks(items, block_size):
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    return [items[i:i + block_size] for i in range(0, len(items), block_size)]


def filter_available_snapshots(start, end, available_snapshots):
    return [numsp for numsp in available_snapshots if start <= numsp <= end]


def format_limited_ints(values, limit=20):
    values = list(values)
    prefix = ",".join(str(value) for value in values[:limit])
    return prefix + ("..." if len(values) > limit else "")


@contextmanager
def exclusive_file_lock(lock_path):
    if fcntl is None:
        yield
        return

    lock_path = Path(lock_path)
    with lock_path.open('w') as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def read_elinfo_csv(path, expected_columns=ELINFO_COLUMNS):
    path = Path(path)
    report = {
        'comments': [],
        'header': None,
        'rows_by_numsp': {},
        'numsp_sequence': [],
        'duplicate_numsp': [],
        'invalid_rows': [],
    }
    if not path.exists():
        return report

    comment_lines = []
    data_lines = []
    with path.open('r', encoding='utf-8', newline='') as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith('#'):
                comment_lines.append(line)
            elif line.strip():
                data_lines.append((line_number, line))

    report['comments'] = comment_lines
    if not data_lines:
        return report

    reader = csv.reader(line for _, line in data_lines)
    rows = list(reader)
    if not rows:
        return report

    header = rows[0]
    report['header'] = header
    missing_columns = [column for column in expected_columns if column not in header]
    if missing_columns:
        raise ValueError(
            "Existing elinfo file is missing summary columns: "
            + ", ".join(missing_columns)
        )

    numsp_index = header.index(NUMSP_COLUMN)
    column_indices = [header.index(column) for column in expected_columns]

    for (line_number, _), row in zip(data_lines[1:], rows[1:]):
        if len(row) != len(header):
            report['invalid_rows'].append(line_number)
            continue
        try:
            numsp = int(float(row[numsp_index]))
        except (TypeError, ValueError):
            report['invalid_rows'].append(line_number)
            continue
        if numsp in report['rows_by_numsp']:
            report['duplicate_numsp'].append(numsp)
        report['numsp_sequence'].append(numsp)
        report['rows_by_numsp'][numsp] = [row[index] for index in column_indices]

    return report


def write_elinfo_csv_atomic(path, comments, rows_by_numsp, columns=ELINFO_COLUMNS):
    path = Path(path)
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            'w',
            encoding='utf-8',
            newline='',
            dir=path.parent,
            prefix=f".{path.name}.tmp.",
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            for comment in comments:
                handle.write(comment if comment.endswith('\n') else comment + '\n')
            writer = csv.writer(handle, lineterminator='\n')
            writer.writerow(columns)
            for numsp in sorted(rows_by_numsp):
                writer.writerow(rows_by_numsp[numsp])
        os.replace(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass


class ProgressReporter:
    def __init__(self, total, desc, enabled=True, unit='snap', mininterval=0.5):
        self.total = total
        self.desc = desc
        self.enabled = enabled and total > 0
        self.unit = unit
        self.mininterval = mininterval
        self.count = 0
        self._last_print = 0.0
        self._use_tqdm = self.enabled and tqdm is not None and sys.stdout.isatty()
        self._progress = None
        if self._use_tqdm:
            self._progress = tqdm(
                total=total,
                desc=desc,
                unit=unit,
                dynamic_ncols=True,
                mininterval=mininterval,
                file=sys.stdout,
            )
        elif self.enabled:
            print(f"[{self.desc}Progress] 0/{self.total} {self.unit}", flush=True)

    def update(self, n=1):
        if not self.enabled:
            return
        self.count += n
        if self._progress is not None:
            self._progress.update(n)
            return

        now = time.perf_counter()
        if self.count >= self.total or now - self._last_print >= self.mininterval:
            print(f"[{self.desc}Progress] {self.count}/{self.total} {self.unit}", flush=True)
            self._last_print = now

    def close(self):
        if self._progress is not None:
            self._progress.close()


def release_process_memory():
    gc.collect()
    if LIBC is not None:
        try:
            LIBC.malloc_trim(0)
        except Exception:
            pass


def get_rss_mb():
    try:
        with open("/proc/self/statm", "r") as fh:
            rss_pages = int(fh.readline().split()[1])
        page_size = os.sysconf("SC_PAGE_SIZE")
        return rss_pages * page_size / (1024.0 * 1024.0)
    except Exception:
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def init_worker(worker_threads):
    global THREADPOOL_CONTROLLER
    if threadpool_limits is not None and worker_threads is not None and worker_threads > 0:
        THREADPOOL_CONTROLLER = threadpool_limits(limits=worker_threads)


class GetInfo:

    def __init__(self):
        model_info = DataProcessor.parse_model_name_details()
        self.modelname = model_info['modelname']
        self.base_modelname = model_info['base_modelname']
        self.dwarf_name = model_info['dwarf_name']
        self.num = model_info['model_num']
        self.description = model_info['description']
        if self.description:
            print(
                f"[GetInfo] Model: {self.base_modelname}, "
                f"description: {self.description}"
            )
        
        # Set up file paths using Path objects
        self.inifile = DataProcessor.resolve_dwarf_ini_path(self.dwarf_name)
        print(f"[GetInfo] Using initial INI file: {self.inifile}")
        
        # Get dwarf parameters
        self.dw_particles = GalaxySimulation.get_dw_num(self.inifile)
        if self.dw_particles <= 0:
            raise ValueError(f"Invalid number of dwarf particles: {self.dw_particles}")
            
        self.core_radius = core_radius_dict.get(self.dwarf_name)
        if self.core_radius is None:
            raise ValueError(f"Core radius not found for dwarf: {self.dwarf_name}")
            
        # Set output file path
        self.output_csv_file = Path(f'./elinfo_{self.modelname}.csv')
        self.profile = False
        self.profile_memory = False
        self.float_precision = 8
        self.show_progress = True

    def get_info(self, numsp):
        start_time = time.perf_counter()
        start_rss = get_rss_mb() if self.profile_memory else None
        snapshot = prepare_snapshot_context(
            folder_path=folder_path,
            snapshot_num=numsp,
            core_radius=self.core_radius,
            include_mw_gas=True,
            mw_gas_radius=300.0,
            include_dark_matter=True,
            include_star_birth=True,
        )
        after_prepare = time.perf_counter()
        after_prepare_rss = get_rss_mb() if self.profile_memory else None

        simulation = snapshot['simulation']
        df = snapshot['df']
        summary = compute_snapshot_summary(snapshot, numsp)
        after_summary = time.perf_counter()
        after_summary_rss = get_rss_mb() if self.profile_memory else None
        simulation.df = None
        del snapshot
        del df, simulation
        release_process_memory()
        end_time = time.perf_counter()

        if self.profile:
            parts = [
                f"pid={os.getpid()}",
                f"numsp={numsp}",
                f"prepare={after_prepare - start_time:.2f}s",
                f"summary={after_summary - after_prepare:.2f}s",
                f"cleanup={end_time - after_summary:.2f}s",
                f"total={end_time - start_time:.2f}s",
            ]
            if self.profile_memory:
                parts.extend([
                    f"rss_start={start_rss:.0f}MB",
                    f"rss_prepare={after_prepare_rss:.0f}MB",
                    f"rss_summary={after_summary_rss:.0f}MB",
                    f"rss_end={get_rss_mb():.0f}MB",
                ])
            print("[GetInfoProfile] " + " ".join(parts), flush=True)
        
        return [summary['star_mass'], summary['star_half_mass'], 
                summary['hotgas_mass'], summary['hotgas_half_mass'], 
                summary['coldgas_mass'], summary['coldgas_half_mass'], 
                summary['mw_mass_r'], summary['eps'], summary['pa'], summary['gas_density'],
                summary['vr'], summary['vtheta'], summary['vphi'], summary['distance_gal'], summary['rhalf'], 
                summary['rhalf_circularized'], summary['rhalf_circular'],
                summary['shape_center_x_kpc'], summary['shape_center_y_kpc'],
                summary['distance'], summary['age'], summary['sigma'], summary['sigma_re_circular'],
                summary['sigma_re_elliptical'], summary['sigma_fixed_500pc'], summary['sigma_re_nstar'],
                summary['sigma_gradient_kms_per_kpc'], summary['tsigma'],
                summary['tsigma_xyz'],
                summary['sigma_x'], summary['sigma_y'], summary['sigma_z'],
                summary['pmra'], summary['pmdec'], summary['cold_gas_center_ra'], summary['cold_gas_center_dec'],
                summary['numsp']]

    def get_info_block(self, numsp_block):
        return [self.get_info(numsp) for numsp in numsp_block]

    def merge_rows_sorted(self, rows, comments):
        formatted_rows = []
        for row in rows:
            numsp = int(row[NUMSP_INDEX])
            formatted_rows.append((
                numsp,
                DataProcessor.format_row_for_csv(row, float_precision=self.float_precision),
            ))

        lock_path = self.output_csv_file.with_name(f".{self.output_csv_file.name}.lock")
        with exclusive_file_lock(lock_path):
            report = read_elinfo_csv(self.output_csv_file)
            output_comments = report['comments'] or comments
            rows_by_numsp = report['rows_by_numsp']
            for numsp, formatted_row in formatted_rows:
                rows_by_numsp[numsp] = formatted_row
            write_elinfo_csv_atomic(self.output_csv_file, output_comments, rows_by_numsp)
        return report

    def warn_about_csv_report(self, report):
        if report['duplicate_numsp']:
            print(
                "GetInfo: removed duplicate numsp row(s) while rewriting CSV: "
                + format_limited_ints(sorted(set(report['duplicate_numsp']))),
                flush=True,
            )
        if report['invalid_rows']:
            print(
                "GetInfo: ignored malformed CSV row(s) while rewriting CSV: line "
                + format_limited_ints(report['invalid_rows']),
                flush=True,
            )

    def process_and_stream(self, pool, numsp_list, comments, chunksize=1, dispatch='blocks', block_size=4):
        total = len(numsp_list)
        progress = ProgressReporter(total, 'GetInfo', enabled=self.show_progress, unit='snap', mininterval=0.5)

        try:
            if dispatch == 'blocks':
                blocks = build_consecutive_blocks(numsp_list, block_size)
                result_iter = pool.imap_unordered(self.get_info_block, blocks, chunksize=1)
                for block_results in result_iter:
                    report = self.merge_rows_sorted(block_results, comments=comments)
                    self.warn_about_csv_report(report)
                    if progress is not None:
                        progress.update(len(block_results))
                    del block_results
                    release_process_memory()
            else:
                result_iter = pool.imap_unordered(self.get_info, numsp_list, chunksize=chunksize)
                for result in result_iter:
                    report = self.merge_rows_sorted([result], comments=comments)
                    self.warn_about_csv_report(report)
                    if progress is not None:
                        progress.update(1)
                    del result
                    release_process_memory()
        finally:
            progress.close()

    def main(self):
        os.environ.setdefault("DSPH_RUN_ID", str(os.getpid()))
        orbit_data = DataProcessor.read_orbit_parameters(orbit_file)
        if orbit_data is None:
            print("Failed to read orbit data.")
            sys.exit(1)

        fd_value = DataProcessor.read_feedback_value(gizmo_config)
        if fd_value is None:
            print("Failed to extract timescale value.")
            sys.exit(1)

        available_snapshots = DataProcessor.list_snapshot_numbers(folder_path)
        if not available_snapshots:
            print(f"No snapshots found under {folder_path}")
            sys.exit(1)

        parser = argparse.ArgumentParser(description='Process some integers.')
        max_snapshot_num = available_snapshots[-1]
        default_range = f'0,{max_snapshot_num}'
        parser.add_argument('--range', type=str, default=default_range, 
                            help='the range of snapshots, e.g., "1,20"')
        parser.add_argument(
            '--output-csv',
            type=Path,
            default=self.output_csv_file,
            help='output elinfo path; defaults to elinfo_<model>.csv',
        )
        parser.add_argument(
            '--rebuild',
            action='store_true',
            help='replace the selected output CSV before processing',
        )
        parser.add_argument('--processes', type=int, default=1, 
                            help='number of processes to use')
        parser.add_argument('--chunksize', type=int, default=1,
                            help='chunksize for pool.imap_unordered')
        parser.add_argument('--dispatch', choices=('blocks', 'dynamic'), default='blocks',
                            help='task dispatch strategy: contiguous snapshot blocks or per-snapshot dynamic scheduling')
        parser.add_argument('--block-size', type=int, default=4,
                            help='number of consecutive snapshots per worker task when --dispatch=blocks')
        parser.add_argument('--maxtasksperchild', type=int, default=4,
                            help='restart each worker after this many pool tasks to reduce memory growth')
        parser.add_argument('--worker-threads', type=int, default=1,
                            help='limit native BLAS/OpenMP threads inside each worker')
        parser.add_argument('--profile', action='store_true',
                            help='print per-snapshot timing information from each worker')
        parser.add_argument('--profile-memory', action='store_true',
                            help='include per-snapshot RSS information in profile logs')
        parser.add_argument('--float-precision', type=int, default=8,
                            help='significant digits for floating-point values written to csv')
        parser.add_argument('--no-progress', action='store_true',
                            help='disable tqdm progress bar output')
        args = parser.parse_args()
        self.profile = args.profile or args.profile_memory
        self.profile_memory = args.profile_memory
        self.float_precision = args.float_precision
        self.show_progress = not args.no_progress
        self.output_csv_file = args.output_csv.resolve()
        if args.rebuild and self.output_csv_file.exists():
            self.output_csv_file.unlink()

        start, end = map(int, args.range.split(','))
        requested_numsp = filter_available_snapshots(start, end, available_snapshots)
        missing_numsp = [numsp for numsp in range(start, end + 1) if numsp not in set(requested_numsp)]
        numsp_list = list(requested_numsp)
        if not requested_numsp:
            print(f"No snapshot files found within requested range {start},{end} under {folder_path}")
            sys.exit(1)
        if missing_numsp:
            print(
                f"GetInfo: skipping {len(missing_numsp)} missing snapshot ids in requested range: "
                + ",".join(str(numsp) for numsp in missing_numsp[:10])
                + ("..." if len(missing_numsp) > 10 else "")
            )

        output_comments = [
            f"# Orbit Data: {','.join(map(str, orbit_data))}\n",
            f"# Timescale_fd: {fd_value}\n",
            "# rhalf: projected elliptical half-light semi-major axis of old stars\n",
            "# rhalf_circularized: equal-area value rhalf*sqrt(1-eps)\n",
            "# rhalf_circular: radius directly enclosing half the old-star light in a circular aperture\n",
            "# sigma and sigma_re_circular: newly formed-star LOS dispersion after removing a planar velocity gradient within circular R<rhalf_circular\n",
            "# sigma_re_elliptical: same tracer and detrending within the old-star half-light ellipse\n",
            "# sigma_fixed_500pc: same tracer and detrending within circular R<0.5 kpc\n",
            "# sigma_re_nstar: newly formed stars in the circular R<rhalf_circular aperture\n",
            "# sigma_gradient_kms_per_kpc: fitted planar LOS velocity-gradient amplitude in that aperture\n",
        ]

        # 断点续写：读取已完成的numsp。读取时会清理历史中断留下的乱序、重复或坏行。
        finished_numsp = set()
        try:
            cleanup_report = self.merge_rows_sorted([], comments=output_comments)
            self.warn_about_csv_report(cleanup_report)
            existing_report = read_elinfo_csv(self.output_csv_file)
            finished_numsp = set(existing_report['rows_by_numsp'])
        except ValueError as e:
            print(str(e))
            print("Please rebuild elinfo from scratch before resuming with the refactored GetInfo.py.")
            sys.exit(1)
        except Exception as e:
            print(f"Warning: Failed to read existing csv: {e}")

        # 只处理未完成的快照
        numsp_list = [n for n in numsp_list if n not in finished_numsp]
        completed_in_range = len(requested_numsp) - len(numsp_list)
        print(
            f"GetInfo resume: total={len(requested_numsp)}, "
            f"completed={completed_in_range}, pending={len(numsp_list)}"
        )

        if numsp_list:
            with get_context("fork").Pool(
                processes=args.processes,
                maxtasksperchild=args.maxtasksperchild,
                initializer=init_worker,
                initargs=(args.worker_threads,),
            ) as pool:
                self.process_and_stream(
                    pool,
                    numsp_list,
                    comments=output_comments,
                    chunksize=args.chunksize,
                    dispatch=args.dispatch,
                    block_size=args.block_size,
                )
            print(f"\nData saved to {self.output_csv_file}")
        else:
            print("All snapshots in the specified range have already been processed.")

        self.sort_csv_by_numsp(comments=output_comments)
        if not self.validate_output(requested_numsp):
            sys.exit(2)

    def sort_csv_by_numsp(self, comments=None):
        report = self.merge_rows_sorted([], comments=comments or [])
        self.warn_about_csv_report(report)

    def validate_output(self, requested_numsp):
        report = read_elinfo_csv(self.output_csv_file)
        expected_numsp = set(requested_numsp)
        found_numsp = set(report['rows_by_numsp'])
        missing_numsp = sorted(expected_numsp - found_numsp)
        sequence = report['numsp_sequence']
        out_of_order = any(left > right for left, right in zip(sequence, sequence[1:]))

        if report['invalid_rows']:
            print(
                "GetInfo validation failed: malformed CSV row(s): line "
                + format_limited_ints(report['invalid_rows'])
            )
            return False
        if report['duplicate_numsp']:
            print(
                "GetInfo validation failed: duplicate numsp row(s): "
                + format_limited_ints(sorted(set(report['duplicate_numsp'])))
            )
            return False
        if out_of_order:
            print("GetInfo validation failed: output CSV is not sorted by numsp.")
            return False
        if missing_numsp:
            print(
                f"GetInfo validation failed: missing {len(missing_numsp)} requested snapshot row(s): "
                + format_limited_ints(missing_numsp)
            )
            return False

        print(
            f"GetInfo validation: complete and sorted "
            f"({len(expected_numsp)} requested snapshot rows present)."
        )
        return True

if __name__ == "__main__":
    processor = GetInfo()
    processor.main()

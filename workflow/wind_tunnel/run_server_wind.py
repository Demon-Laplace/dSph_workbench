"""Stage a frozen wind-only case, then run it in one astro_thin allocation.

Staging copies existing inputs only. Solver compilation, integration and
scientific measurement run inside Slurm; no local simulation is launched.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import time

from run_local_wind import replace_parameters, THREAD_ENV


def utc():
    return datetime.now(timezone.utc).isoformat()


def write(path, value):
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n', encoding='utf8')
    temporary.replace(path)


def stage(args):
    project = args.project.resolve()
    local = json.loads(args.source_runtime.read_text(encoding='utf8'))
    source = project/local['source_case']
    origin = local['wind_origin_time_code']
    parameters = (project/local['case']/'GZMB_LONG.PARAM').read_text(encoding='utf8')
    parameters = replace_parameters(parameters, {
        'InitCondFile': './initial_000', 'OutputDir': './output/',
        'TimeBegin': format(origin, '.17g'), 'TimeMax': format(local['time_max_code'], '.17g'),
        'TimeLimitCPU': 410400, 'CpuTimeBetRestartFile': 1800,
        'NumFilesWrittenInParallel': 4, 'MaxMemSize': 2500,
        'PartAllocFactor': 2.5, 'BufferSize': 95,
        'TreeDomainUpdateFrequency': '.002'})
    config = {
        'model': local['model'], 'case': local['model'],
        'solver': 'solver/GIZMO_wind_astrothin', 'parameters': 'GZMB_LONG.PARAM',
        'initial': local['initial'], 'forcing': local['forcing'],
        'frozen_solver_archive': local['solver_archive'],
        'wind_origin_time_code': origin, 'time_max_code': local['time_max_code'],
        'target_elapsed_Myr': local['target_elapsed_Myr'], 'fresh_from_initial': True,
        'predecessor': local['case'], 'predecessor_resumed': False,
        'mpi_processes': 24, 'threads_per_process': 1, 'nodes': 1,
        'cluster': 'astro_thin', 'gpu_enabled': False,
        'galactic_gravity': False, 'isolation_comparison': False,
        'common_core_translation_and_boost': True,
        'technical_overrides': {'BufferSize': 95, 'MaxMemSize': 2500,
                                'PartAllocFactor': 2.5, 'NumFilesWrittenInParallel': 4,
                                'TreeDomainUpdateFrequency': .002},
        'created_at_utc': utc(), 'result_stem': 'wind_only_astrothin'}
    temporary = args.output.parent/(args.output.stem+'_stage')
    temporary.mkdir(exist_ok=False)
    try:
        (temporary/'GZMB_LONG.PARAM').write_text(parameters, encoding='utf8')
        write(temporary/'run.json', config)
        tools = Path(__file__).resolve().parent
        with tarfile.open(args.output, 'w') as archive:
            archive.add(temporary/'run.json', arcname='run.json')
            archive.add(temporary/'GZMB_LONG.PARAM', arcname=local['model']+'/GZMB_LONG.PARAM')
            for name in ('TREECOOL', 'ewald_spc_table_64_dbl.dat'):
                archive.add(source/name, arcname=local['model']+'/'+name)
            archive.add(project/local['initial'], arcname=local['model']+'/initial_000.hdf5')
            archive.add(project/local['forcing'], arcname=local['model']+'/wind_tunnel_orbit.dat')
            for name in ('wind_time_origin.dat', 'output_times.txt'):
                archive.add(project/local['case']/name, arcname=local['model']+'/'+name)
            archive.add(project/local['solver_archive'], arcname='solver_source.tar.gz')
            for name in ('run_server_wind.py', 'run_local_wind.py'):
                archive.add(tools/name, arcname='tools/'+name)
            archive.add(tools/'wind_only_astrothin.slurm', arcname='wind_only_astrothin.slurm')
            archive.add(project/'Code/analysis/model_comparisons/wind_tunnel_response.py',
                        arcname='tools/wind_tunnel_response.py')
        write(args.manifest, config | {'bundle': args.output.name, 'bundle_bytes': args.output.stat().st_size})
    finally:
        shutil.rmtree(temporary)
    print(json.dumps({'bundle': str(args.output), 'bytes': args.output.stat().st_size}))


def last_sync(path):
    if not path.exists():
        return None
    with path.open('rb') as stream:
        stream.seek(max(0, path.stat().st_size-65536))
        text = stream.read().decode('ascii', errors='replace')
    rows = re.findall(r'Sync-Point\s+(\d+),\s+Time:\s*([0-9.eE+-]+)', text)
    return (int(rows[-1][0]), float(rows[-1][1])) if rows else None


def rank_affinity(solver, case):
    found = []
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():
            continue
        try:
            executable = (entry/'cmdline').read_bytes().split(b'\0')[0].decode()
            if executable != str(solver) or (entry/'cwd').resolve() != case:
                continue
            text = (entry/'status').read_text()
            found.append({'pid': int(entry.name),
                          'cpus_allowed': re.search(r'Cpus_allowed_list:\s*(\S+)', text)[1]})
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return found


def build(root):
    with tarfile.open(root/'solver_source.tar.gz') as archive:
        for member in archive.getmembers():
            if not (root/member.name).resolve().is_relative_to(root):
                raise ValueError('Solver archive contains an external path')
        archive.extractall(root)
    wrapper = shutil.which('h5cc')
    if not wrapper:
        raise RuntimeError('The registered HDF5 module does not supply h5cc')
    hdf = Path(wrapper).resolve().parent.parent
    command = ['make', '-B', '-j4', 'SYSTYPE=Meso', 'CC=mpicc', 'FC=mpicc',
               'EXEC=GIZMO_wind_astrothin',
               'OPTIMIZE=-O3 -march=native -std=gnu11 -fcommon -ffp-contract=off',
               'HDF5INCL=-I'+str(hdf/'include')+' -DH5_USE_16_API',
               'HDF5LIB=-L'+str(hdf/'lib')+' -Wl,-rpath,'+str(hdf/'lib')+' -lhdf5 -lz',
               'MPICHLIB=', 'GSL_LIBS=', 'FFTW_LIBS=', 'FFTW_LIBNAMES=', 'GRACKLELIBS=']
    with (root/'build.log').open('w') as output:
        subprocess.run(command, cwd=root/'solver', stdout=output, stderr=subprocess.STDOUT, check=True)


def run(args):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Compile and execute only inside Slurm')
    root = args.config.resolve().parent
    config = json.loads(args.config.read_text(encoding='utf8'))
    ranks = int(os.environ['SLURM_NTASKS'])
    if ranks != config['mpi_processes'] or int(os.environ.get('SLURM_JOB_NUM_NODES', '1')) != 1:
        raise ValueError('This shared-memory case requires its declared single-node MPI layout')
    os.environ.update(THREAD_ENV)
    case = (root/config['case']).resolve()
    (case/'output').mkdir(exist_ok=True)
    solver = (root/config['solver']).resolve()
    if args.mode == 2:
        if list((case/'output').glob('snapshot_*.hdf5')):
            raise ValueError('Fresh launch would overwrite evolved output; use native restart')
        build(root)
    elif not solver.exists():
        raise FileNotFoundError('Resume with the existing compatible frozen executable')
    baseline = last_sync(case/'output/info.txt')
    initial_time = baseline[1] if baseline else config['wind_origin_time_code']
    command = ['mpirun', '-np', str(ranks), str(solver), config['parameters'], str(args.mode)]
    status = config | {'status': 'running', 'job_id': os.environ['SLURM_JOB_ID'],
                       'node': os.uname().nodename, 'started_at_utc': utc(),
                       'restart_flag': args.mode, 'command': command}
    write(root/'status.json', status)
    started = time.monotonic()
    with (case/'solver.stdout.txt').open('a') as output:
        process = subprocess.Popen(command, cwd=case, stdout=output, stderr=subprocess.STDOUT)
        while process.poll() is None:
            state = last_sync(case/'output/info.txt')
            if time.monotonic()-started >= 140 and state and state[1]>initial_time+1.e-8:
                write(root/'startup.json', status | {'status': 'startup_passed',
                      'detected_at_utc': utc(), 'running_seconds': time.monotonic()-started,
                      'latest_sync': state, 'rank_affinity': rank_affinity(solver, case)})
                break
            time.sleep(15)
        code = process.wait()
    state = last_sync(case/'output/info.txt')
    reached = bool(state and state[1] >= config['time_max_code']-1.e-6)
    status.update(status='completed' if code == 0 and reached else ('checkpointed' if code == 0 else 'failed'),
                  exit_code=code, latest_sync=state, finished_at_utc=utc())
    write(root/'status.json', status)
    if code:
        raise SystemExit(code)
    measure = root/'tools/wind_tunnel_response.py'
    subprocess.run([sys.executable, str(measure), 'measure', '--kind', 'wind',
                    '--model', config['model'], '--directory', str(case/'output'),
                    '--initial', str(case/'initial_000.hdf5'), '--forcing', str(case/'wind_tunnel_orbit.dat'),
                    '--output', str(case/'response.csv'), '--append'], check=True)
    subprocess.run([sys.executable, str(measure), 'summary', '--wind', str(case/'response.csv'),
                    '--wind-only', '--duration-myr', str(config['target_elapsed_Myr']),
                    '--output', str(case/'summary.json')], check=True)
    subprocess.run([sys.executable, str(measure), 'plot', '--wind-series', config['model'],
                    str(case/'response.csv'), '--wind-only', '--overview',
                    '--output', str(case/(config['model']+'_wind_only_timeseries.pdf'))], check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    s = commands.add_parser('stage')
    s.add_argument('--project', type=Path, required=True)
    s.add_argument('--source-runtime', type=Path, required=True)
    s.add_argument('--output', type=Path, required=True)
    s.add_argument('--manifest', type=Path, required=True)
    r = commands.add_parser('run')
    r.add_argument('--config', type=Path, required=True)
    r.add_argument('--mode', type=int, choices=(1, 2), default=2)
    args = parser.parse_args()
    stage(args) if args.command == 'stage' else run(args)

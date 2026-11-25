"""Continue a local wind-only case, publishing closed outputs to the project disk.

The first segment can redistribute a saved HDF5 state onto a new MPI layout.
Later segments use native checkpoints. Run this entry point inside Linux; the
workspace must be on local Linux storage, separate from the project archive.
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
import traceback


TIME_MYR = 3.08568e16 / 3.15576e13
THREAD_ENV = {name: "1" for name in
              ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def replace_parameters(text, changes):
    for key, value in changes.items():
        pattern = rf"(?m)^{re.escape(key)}\s+[^\n]+$"
        text, count = re.subn(pattern, lambda _: f"{key:<34}{value}", text)
        if count != 1:
            raise ValueError(f"Expected one parameter {key}, got {count}")
    return text


def snapshot_time(path):
    import h5py
    with h5py.File(path) as handle:
        return float(handle["Header"].attrs["Time"])


def archive_solver(solver, destination):
    with tarfile.open(destination, "w:gz") as archive:
        def keep(member):
            return None if member.name.endswith((".o", ".oo")) else member
        archive.add(solver, arcname="solver", filter=keep)


def stage(args):
    project = args.project.resolve()
    case = (project / args.case).resolve()
    source = (project / (args.source_case or args.case)).resolve()
    workspace = args.workspace.resolve()
    if str(workspace).startswith("/mnt/") or workspace == Path("/"):
        raise ValueError("The active workspace must be on local Linux storage")
    if (workspace / "case").exists():
        raise FileExistsError("An existing case must be resumed, not restaged")
    if args.from_initial and source == case:
        raise ValueError("A fresh round needs a distinct archive case")
    if source != case:
        case.mkdir(parents=True, exist_ok=False)
        (case / "output").mkdir()
    workspace.mkdir(parents=True, exist_ok=True)
    run = workspace / "case"
    run.mkdir()
    (run / "output").mkdir()
    (run / "restart_input").mkdir()
    initial = (project / args.initial).resolve() if args.initial else source / "initial_000.hdf5"
    forcing = (project / args.forcing).resolve() if args.forcing else source / "wind_tunnel_orbit.dat"
    latest = (initial if args.from_initial else
              sorted((source / "output").glob("snapshot_[0-9][0-9][0-9].hdf5"))[-1])
    start = snapshot_time(initial)
    current = snapshot_time(latest)
    shutil.copyfile(latest, run / "restart_input" / latest.name)
    for name in ("TREECOOL", "ewald_spc_table_64_dbl.dat"):
        shutil.copyfile(source / name, run / name)
    shutil.copyfile(forcing, run / "wind_tunnel_orbit.dat")
    (run / "wind_time_origin.dat").write_text(f"{start:.17g}\n")
    target_time = start + args.target_myr / TIME_MYR
    changes = {
        "InitCondFile": f"./restart_input/{latest.stem}",
        "RestartFile": "long_restart",
        "OutputListOn": 0,
        "TimeBegin": f"{current:.17g}",
        "TimeMax": f"{target_time:.17g}",
        "TimeOfFirstSnapshot": f"{start:.17g}",
        "TimeBetSnapshot": ".005",
        "TimeLimitCPU": args.first_segment_seconds,
        "CpuTimeBetRestartFile": 1800,
    }
    if args.domain_frequency is not None:
        changes["TreeDomainUpdateFrequency"] = args.domain_frequency
    if args.buffer_mib is not None:
        changes["BufferSize"] = args.buffer_mib
    if args.from_initial:
        # Keep the first pilot diagnostic and the regular cadence, while the
        # integer clock always spans the full planned integration interval.
        count = int((target_time-start)/.005)
        outputs = [start+.001] + [start+i*.005 for i in range(1, count+1)]
        outputs = [value for value in outputs if value < target_time-1.e-10] + [target_time]
        if len(outputs) > 1200:
            raise ValueError("Output list exceeds this solver's MAXLEN_OUTPUTLIST")
        output_text = "".join(f"{value:.17g}\n" for value in outputs)
        (run / "output_times.txt").write_text(output_text)
        (case / "output_times.txt").write_text(output_text)
        changes["OutputListOn"] = 1
    params = replace_parameters((source / "GZMB.PARAM").read_text(), changes)
    (run / "GZMB_LONG.PARAM").write_text(params)
    shutil.copyfile(run / "GZMB_LONG.PARAM", case / "GZMB_LONG.PARAM")
    shutil.copyfile(run / "wind_time_origin.dat", case / "wind_time_origin.dat")
    model_dir = case.parent
    support = model_dir / "_support"
    stem = args.result_stem
    short_series = support / "wind_only_local_response.csv"
    series = support / f"{stem}_response.csv"
    if not args.from_initial and source == case and not series.exists() and short_series.exists():
        shutil.copyfile(short_series, series)
    config = {
        "model": args.model, "case": str(case.relative_to(project)),
        "source_case": str(source.relative_to(project)),
        "initial": str(initial.relative_to(project)), "forcing": str(forcing.relative_to(project)),
        "fresh_from_initial": args.from_initial,
        "project": str(project), "workspace": str(workspace),
        "target_elapsed_Myr": args.target_myr, "wind_origin_time_code": start,
        "time_begin_code": current, "time_max_code": target_time,
        "resumed_snapshot": latest.name, "resumed_elapsed_Myr": (current-start)*TIME_MYR,
        "next_restart_flag": 2, "mpi_processes": len(args.cpus), "cpu_affinity": args.cpus,
        "nice_priority": 10, "threads_per_process": 1,
        "segment_wall_limit_seconds": args.segment_seconds,
        "first_segment_wall_limit_seconds": args.first_segment_seconds,
        "snapshot_interval_Myr": .005*TIME_MYR,
        "result_stem": stem, "figure": args.figure or f"{args.model}_wind_only_timeseries.pdf",
        "technical_overrides": {key: changes[key] for key in ("TreeDomainUpdateFrequency", "BufferSize") if key in changes},
        "solver": str(workspace / "solver" / "GIZMO_wind"),
        "solver_archive": str((support / f"{stem}_solver.tar.gz").relative_to(project)),
        "galactic_gravity": False, "isolation_comparison": False,
        "common_core_translation_and_boost": True,
        "created_at_utc": utc(),
    }
    cfg = support / f"{stem}_runtime.json"
    write_json(cfg, config)
    # Freeze the derived solver, excluding disposable build objects.
    archive_solver(workspace / "solver", support / f"{stem}_solver.tar.gz")
    env = os.environ | THREAD_ENV
    log = case / "long_supervisor.txt"
    with log.open("a") as out:
        process = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "supervise", "--config", str(cfg)],
            stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT,
            start_new_session=True, env=env)
    print(json.dumps({"supervisor_pid": process.pid, "config": str(cfg),
                      "continued_elapsed_Myr": config["resumed_elapsed_Myr"]}))


def mirror_snapshots(run, case, completed):
    for number in sorted(completed):
        source = run / "output" / f"snapshot_{number:03d}.hdf5"
        destination = case / "output" / source.name
        if not source.exists():
            continue
        if destination.exists():
            raise FileExistsError(f"Refusing to replace archived snapshot {destination}")
        partial = destination.with_suffix(".partial.hdf5")
        shutil.copyfile(source, partial)
        partial.replace(destination)
        source.unlink()
        print(f"{utc()} published {destination.name}", flush=True)
    completed.clear()


def publish_logs(run, case, *, checkpoints=False):
    for path in (run / "output").iterdir():
        if path.is_file() and not path.name.startswith("snapshot_") and path.name not in ("stop", "cont"):
            temporary = case / "output" / ("long_" + path.name + ".partial")
            shutil.copyfile(path, temporary)
            temporary.replace(temporary.with_name(temporary.name.removesuffix(".partial")))
    if checkpoints:
        destination = case / "output" / "restartfiles"
        destination.mkdir(exist_ok=True)
        for path in (run / "output" / "restartfiles").glob("long_restart.*"):
            if path.name.endswith(".bak"):
                continue
            temporary = destination / (path.name + ".partial")
            shutil.copyfile(path, temporary)
            temporary.replace(destination / path.name)
    for path in run.glob("GZMB_LONG.PARAM*"):
        shutil.copyfile(path, case / path.name)


def live_clock(run, default):
    path = run / "output" / "info.txt"
    if not path.exists():
        return default
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size-4096))
        tail = handle.read().decode(errors="replace")
    matches = re.findall(r"Time:\s*([\deE.+-]+)", tail)
    return float(matches[-1]) if matches else default


def diagnostics(config, status):
    project = Path(config["project"])
    case = project / config["case"]
    support = case.parent / "_support"
    stem = config["result_stem"]
    processor = project / "Code/analysis/model_comparisons/wind_tunnel_response.py"
    series = support / f"{stem}_response.csv"
    initial = project / config.get("initial", str(Path(config["case"]) / "initial_000.hdf5"))
    forcing = project / config.get("forcing", str(Path(config["case"]) / "wind_tunnel_orbit.dat"))
    common = ["nice", "-n", "10", "taskset", "-c", ",".join(map(str, config["cpu_affinity"])),
              sys.executable, str(processor)]
    commands = [
        ["measure", "--kind", "wind", "--model", config["model"],
         "--directory", str(case / "output"), "--initial", str(initial),
         "--forcing", str(forcing), "--output", str(series), "--append"],
        ["summary", "--wind", str(series), "--output", str(support / f"{stem}_summary.json"),
         "--duration-myr", str(config["target_elapsed_Myr"]), "--wind-only"],
    ]
    temporary_pdf = case.parent / (Path(config["figure"]).stem + ".partial.pdf")
    commands.append(["plot", "--wind-series", config["model"], str(series),
                     "--wind-only", "--overview", "--output", str(temporary_pdf)])
    env = os.environ | THREAD_ENV | {"MPLBACKEND": "Agg", "MPLCONFIGDIR": str(Path(config["workspace"]) / "mpl")}
    for command in commands:
        subprocess.run(common + command, check=True, env=env)
    temporary_pdf.replace(case.parent / config["figure"])
    status["diagnostics_updated_at_utc"] = utc()
    status.pop("diagnostics_error", None)


def supervise(args):
    config = json.loads(args.config.read_text())
    project = Path(config["project"])
    case = project / config["case"]
    workspace = Path(config["workspace"])
    run = workspace / "case"
    support = case.parent / "_support"
    status_path = support / f"{config['result_stem']}_status.json"
    stop_path = support / f"{config['result_stem']}.stop"
    status = {"status": "starting", "supervisor_pid": os.getpid(),
              "started_at_utc": utc(), "target_elapsed_Myr": config["target_elapsed_Myr"],
              "elapsed_Myr": config["resumed_elapsed_Myr"], "segment": 0,
              "mpi_processes": config["mpi_processes"], "cpu_affinity": config["cpu_affinity"],
              "galactic_gravity": False, "isolation_comparison": False}
    write_json(status_path, status)
    stdout_path = case / "long_solver.stdout.txt"
    try:
        while not stop_path.exists():
            status["segment"] += 1
            status["status"] = "running"
            command = ["nice", "-n", str(config["nice_priority"]), "mpirun", "-np",
                       str(config["mpi_processes"]), "-bind-to",
                       "user:" + ",".join(map(str, config["cpu_affinity"])),
                       config["solver"], "GZMB_LONG.PARAM", str(config["next_restart_flag"])]
            completed, writing, buffer = set(), None, ""
            final_reached = False
            segment_started = time.monotonic()
            with stdout_path.open("ab") as out:
                offset = out.tell()
                process = subprocess.Popen(command, cwd=run, stdout=out, stderr=subprocess.STDOUT,
                                           env=os.environ | THREAD_ENV)
                status.update(mpi_launcher_pid=process.pid, mpi_command=command,
                              segment_started_at_utc=utc())
                write_json(status_path, status)
                with stdout_path.open("r", errors="replace") as reader:
                    reader.seek(offset)
                    while True:
                        exit_code = process.poll()
                        buffer += reader.read()
                        lines = buffer.split("\n")
                        buffer = lines.pop()
                        for line in lines:
                            match = re.search(r"writing snapshot file #(\d+)", line)
                            if match:
                                writing = int(match[1])
                            if "done with snapshot." in line and writing is not None:
                                completed.add(writing)
                                status["last_completed_snapshot"] = writing
                                writing = None
                            if "Simulation ends." in line:
                                final_reached = True
                        mirror_snapshots(run, case, completed)
                        fallback = config["wind_origin_time_code"] + status["elapsed_Myr"]/TIME_MYR
                        clock = live_clock(run, fallback)
                        status.update(time_code=clock,
                                      elapsed_Myr=(clock-config["wind_origin_time_code"])*TIME_MYR,
                                      updated_at_utc=utc(),
                                      segment_wall_seconds=time.monotonic()-segment_started)
                        publish_logs(run, case)
                        write_json(status_path, status)
                        if exit_code is not None:
                            break
                        if stop_path.exists():
                            (run / "output" / "stop").touch()
                        time.sleep(10)
            status["last_exit_code"] = process.returncode
            publish_logs(run, case, checkpoints=True)
            if process.returncode:
                raise RuntimeError(f"MPI integration exited {process.returncode}; preserve the local checkpoint")
            config["next_restart_flag"] = 1
            params = run / "GZMB_LONG.PARAM"
            params.write_text(replace_parameters(params.read_text(), {
                "TimeLimitCPU": config["segment_wall_limit_seconds"]}))
            shutil.copyfile(params, case / params.name)
            write_json(args.config, config)
            status["status"] = "updating_diagnostics"
            write_json(status_path, status)
            try:
                diagnostics(config, status)
            except Exception as error:
                status["diagnostics_error"] = str(error)
                traceback.print_exc()
            if final_reached:
                status.update(status="completed", elapsed_Myr=config["target_elapsed_Myr"],
                              finished_at_utc=utc())
                # All completed snapshots and the latest native checkpoint are
                # already published. Only this task's local workspace is removed.
                if workspace.parent != Path("/opt") or not workspace.name.startswith("f8002w-"):
                    raise ValueError("Unexpected cleanup root; leave the local workspace for manual review")
                shutil.rmtree(workspace)
                status["local_workspace_removed"] = True
                break
            if stop_path.exists():
                status["status"] = "stopped"
                break
            if not (run / "output" / "cont").exists():
                raise RuntimeError("Normal exit without final time, stop request or wall-time continuation marker")
            (run / "output" / "cont").unlink()
            status["status"] = "restarting_native_checkpoint"
            write_json(status_path, status)
        else:
            status["status"] = "stopped"
    except Exception as error:
        status.update(status="failed", error=str(error), failed_at_utc=utc())
        traceback.print_exc()
    finally:
        status["updated_at_utc"] = utc()
        write_json(status_path, status)
        print(json.dumps(status, ensure_ascii=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    start = commands.add_parser("start")
    start.add_argument("--project", type=Path, required=True)
    start.add_argument("--case", type=Path, required=True, help="Relative to the project")
    start.add_argument("--source-case", type=Path, help="Existing case supplying validated parameters and inputs")
    start.add_argument("--initial", type=Path, help="Shared initial snapshot, relative to the project")
    start.add_argument("--forcing", type=Path, help="Shared forcing table, relative to the project")
    start.add_argument("--from-initial", action="store_true", help="Start a fresh round in a distinct archive case")
    start.add_argument("--domain-frequency", type=float, help="Override only the domain decomposition threshold")
    start.add_argument("--buffer-mib", type=int, help="Override the communication buffer per MPI rank")
    start.add_argument("--model", required=True)
    start.add_argument("--workspace", type=Path, required=True, help="Prepared local solver workspace")
    start.add_argument("--target-myr", type=float, required=True)
    start.add_argument("--cpus", type=int, nargs="+", required=True)
    start.add_argument("--segment-seconds", type=int, default=7200)
    start.add_argument("--first-segment-seconds", type=int, default=360,
                       help="First checkpoint interval for verifying the new MPI layout")
    start.add_argument("--result-stem", default="wind_only_long")
    start.add_argument("--figure", help="Output PDF filename under this model")
    resume = commands.add_parser("supervise")
    resume.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "start":
        stage(args)
    else:
        supervise(args)


if __name__ == "__main__":
    main()

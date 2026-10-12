#!/usr/bin/env python3
"""Inventory, validate, and prepare isolated analysis runs; never starts simulations."""

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import workflow_runtime as runtime
import workflow_management as management
import workflow_models as models
import workflow_release as release
import workflow_github as github

ROOT = Path(__file__).resolve().parent
MODEL = re.compile(r"(?:Fornax|MW)\d+(?:_.+)?$")
SKIP = {".git", ".aws", ".codex", ".agents", "__pycache__", ".ipynb_checkpoints"}


def files_under(root, exclude=()):
    """Walk regular files without following directory or file symlinks."""
    for base, dirs, files in os.walk(root):
        dirs[:] = sorted(d for d in dirs if d not in SKIP and d not in exclude
                         and not (Path(base) / d).is_symlink())
        for name in sorted(files):
            path = Path(base) / name
            if path.is_file() and not path.is_symlink():
                yield path


def model_rows(root):
    for directory in sorted(root.iterdir()):
        if not directory.is_dir() or directory.is_symlink() or not MODEL.fullmatch(directory.name):
            continue
        files = list(files_under(directory))
        snapshots = [p for p in files if re.fullmatch(r"snapshot_\d+(?:\.\d+)?\.hdf5", p.name)]
        yield {
            "model": directory.name,
            "files": len(files),
            "GiB": round(sum(p.stat().st_size for p in files) / 2**30, 4),
            "snapshot_files": len(snapshots),
            "snapshot_locations": ";".join(sorted({str(p.parent.relative_to(root)) for p in snapshots})),
            "videos": sum(p.suffix == ".mp4" for p in files),
            "images": sum(p.suffix == ".png" for p in files),
            "summary_tables": ";".join(str(p.relative_to(root)) for p in files if p.name.startswith("elinfo_") and p.suffix == ".csv"),
            "local_scripts": sum(p.suffix == ".py" for p in files),
            "has_orbit": (directory / "mergegalaxy.par").is_file(),
            "has_gizmo_parameters": (directory / "GZWJL.PARAM").is_file(),
        }


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory(args):
    rows = list(model_rows(args.root))
    if args.model:
        rows = [r for r in rows if r["model"] == args.model]
    if args.csv:
        if not rows:
            raise ValueError("No matching model directories")
        writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    else:
        print(f"{'Model':25} {'GiB':>8} {'Snaps':>6} {'Videos':>6} {'Scripts':>7}  Snapshot locations")
        for row in rows:
            print(f"{row['model']:25} {row['GiB']:8.3f} {row['snapshot_files']:6} "
                  f"{row['videos']:6} {row['local_scripts']:7}  {row['snapshot_locations'] or '-'}")
        print(f"\n{len(rows)} model directories; {sum(r['GiB'] for r in rows):.3f} GiB inside model directories.")
        print("Snapshots count files, including separate copies/parts; root-level shared inputs are excluded.")
    return 0


def audit(args):
    upstream = args.repo / "repo"
    if not upstream.is_dir():
        raise ValueError(f"Missing upstream analysis directory: {upstream}")
    candidates = {}
    for path in files_under(upstream):
        if path.suffix == ".py":
            candidates.setdefault(path.name, []).append(path)
    print("Local script\tComparison with upstream\tUpstream path")
    for path in files_under(args.root, {'workspace', 'tools', 'runs', 'analysis'}):
        if path.suffix != ".py" or args.repo in path.parents or path.name == "workflow.py":
            continue
        # Generated runs and independent checkouts are not legacy source trees.
        if path.relative_to(args.root).parts[0] in {"tools", "runs", "analysis"}:
            continue
        matches = candidates.get(path.name, [])
        status = "local only"
        if len(matches) == 1:
            status = "identical" if digest(path) == digest(matches[0]) else "different; review before replacing"
        elif len(matches) > 1:
            status = "ambiguous basename; review manually"
        print(f"{path.relative_to(args.root)}\t{status}\t{';'.join(str(p.relative_to(args.repo)) for p in matches)}")
    return 0


def git_output(directory, *arguments):
    result = subprocess.run(["git", "-C", str(directory), *arguments], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def doctor(args):
    print("Workspace:", args.root)
    print("Python:", sys.executable)
    print("Analysis checkout:", args.repo)
    print("Analysis commit:", git_output(args.repo, "rev-parse", "HEAD"))
    print("Workspace Git HEAD:", git_output(args.root, "rev-parse", "--verify", "HEAD"))
    missing = []
    for module in ("numpy", "scipy", "pandas", "matplotlib", "astropy", "h5py", "photutils", "sklearn"):
        found = importlib.util.find_spec(module) is not None
        print(f"{module}: {'installed (import not tested)' if found else 'MISSING'}")
        if not found:
            missing.append(module)
    print("agama (IC generation):", "installed" if importlib.util.find_spec("agama") else "MISSING")
    print("Upstream model_tool infers the model prefix from its directory name.")
    print("Use check-ini below for validation by explicit path; direct builds also auto-push to GitHub.")
    print("Configured production code:", runtime.code_dir(args))
    print("Configuration:", args.root / "workflow.json")
    return int(bool(missing) or not (args.repo / "repo" / "GetInfo.py").is_file())


def check_ini(args):
    script = args.repo / "model_tool.py"
    if not script.is_file():
        raise ValueError(f"Missing model tool: {script}")
    paths = []
    for selector in args.models:
        for token in selector.split(","):
            token = token.strip()
            match = re.fullmatch(r"(\d+)(?:-(\d+))?", token)
            if not match:
                raise ValueError(f"Expected model numbers or inclusive ranges, got {token!r}")
            start = int(match[1])
            end = int(match[2] or match[1])
            if start > end or end - start > 10000:
                raise ValueError(f"Invalid or excessively large range: {token}")
            paths.extend(args.repo / "IniFile" / f"IC_Fornax{number}.ini" for number in range(start, end + 1))
    paths = list(dict.fromkeys(paths))
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise ValueError("Missing parameter files: " + ", ".join(missing))
    return subprocess.run([sys.executable, str(script), "check", *map(str, paths)]).returncode


def prepare(args):
    return runtime.prepare(args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="simulation workspace (default: this script's directory)")
    parser.add_argument("--repo", type=Path, help="dSph_workbench checkout (default: ROOT/tools/dSph_workbench)")
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("sync-github", help="validate, commit and push production source to the dedicated GitHub branch")
    p.add_argument("--message", required=True)
    p.set_defaults(func=github.sync)
    for name, function in (("release", release.build), ("deploy", release.deploy),
                           ("activate-release", release.activate)):
        p = commands.add_parser(name, help="create, transfer or activate a checksum-verified Git release")
        p.add_argument("--name", required=True)
        p.set_defaults(func=function)
    p = commands.add_parser("list", help="inventory local model data without reading snapshot contents")
    p.add_argument("--model", help="one exact directory name, e.g. Fornax1071")
    p.add_argument("--csv", action="store_true", help="write CSV to stdout")
    p.set_defaults(func=inventory)
    p = commands.add_parser("audit", help="compare local Python copies with GitHub checkout; no replacements")
    p.set_defaults(func=audit)
    p = commands.add_parser("doctor", help="inspect paths, code revision, and installed dependency availability")
    p.set_defaults(func=doctor)
    p = commands.add_parser("check-ini", help="validate GitHub parameter files without writes or automatic push")
    p.add_argument("models", nargs="+", help='e.g. "4001,7002-7003"')
    p.set_defaults(func=check_ini)
    p = commands.add_parser("prepare", help="prepare a separate analysis directory with pinned code and linked data")
    p.add_argument("model", help="e.g. Fornax1071")
    p.add_argument("--name", required=True, help="unique analysis run name")
    p.add_argument("--snapshots", choices=("model", "output"), help="explicit choice when snapshots exist in both locations")
    p.add_argument("--dry-run", action="store_true", help="validate and display paths without creating files")
    p.set_defaults(func=prepare)
    p = commands.add_parser("new-ini", help="create a draft parameter file without overwriting GitHub or pushing")
    p.add_argument("model")
    p.add_argument("--base", required=True)
    p.add_argument("--set", dest="settings", action="append", default=[])
    p.set_defaults(func=management.new_ini)
    p = commands.add_parser("push-params", help="upload a versioned parameter batch and verify SHA256 on the server")
    p.add_argument("models", nargs="+")
    p.add_argument("--batch", required=True)
    p.set_defaults(func=management.push_params)
    p = commands.add_parser("verify-params", help="server: verify received parameter batch")
    p.add_argument("--batch", required=True)
    p.set_defaults(func=management.verify_params)
    p = commands.add_parser("fetch-model", help="retrieve existing model figures/videos/tables without snapshots")
    p.add_argument("model")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(func=management.fetch_model)
    p = commands.add_parser("new-case", help="prepare a new production case from an explicit template and particle IC")
    p.add_argument("model")
    p.add_argument("--template", required=True)
    p.add_argument("--dwarf-ic", required=True, help="particle IC path, relative to workspace or absolute")
    p.add_argument("--mw-ic", help="override template Milky Way particle IC")
    p.set_defaults(func=models.new_case)
    p = commands.add_parser("simulate", help="preflight or explicitly submit a prepared simulation; refuses existing data")
    p.add_argument("model")
    p.add_argument("--name", required=True)
    p.add_argument("--launcher", choices=("subastrothin.sh", "subastrofat.sh", "submeso1.sh"), default="subastrothin.sh")
    p.add_argument("--submit", action="store_true")
    p.set_defaults(func=models.simulate)
    p = commands.add_parser("collect", help="wait for a remote analysis and automatically retrieve results")
    p.add_argument("--name", required=True)
    p.add_argument("--wait", action="store_true")
    p.add_argument("--timeout", type=int, default=7200, help="maximum wait in seconds")
    p.set_defaults(func=management.collect)
    for name, function, help_text in (
        ("organize", management.organize, "refresh categorized workspace and catalogue"),
        ("pull", management.pull, "download analysis outputs without snapshots"),
        ("status", management.status, "show run state and server queues"),
        ("remote", management.remote, "execute this workflow on the configured server"),
        ("run", runtime.run, "run and validate analysis in a prepared directory"),
        ("submit", runtime.submit, "submit a prepared analysis to Slurm"),
        ("pipeline", management.pipeline, "prepare and schedule an isolated analysis"),
    ):
        p = commands.add_parser(name, help=help_text)
        p.set_defaults(func=function)
        if name == "remote":
            p.add_argument("arguments", nargs=argparse.REMAINDER)
        elif name == "status":
            p.add_argument("--name")
            p.add_argument("--json", action="store_true")
        elif name != "organize":
            p.add_argument("--name", required=True)
        if name in ("pull", "pipeline"):
            p.add_argument("--dry-run", action="store_true")
        if name == "pipeline":
            p.add_argument("model")
            p.add_argument("--snapshots", choices=("model", "output"))
            p.add_argument("--submit", action="store_true", help="submit job; otherwise only validate with sbatch --test-only")
        if name in ("run", "submit", "pipeline"):
            p.add_argument("--range", required=True, help="inclusive START,END")
            p.add_argument("--processes", type=int, default=2)
            p.add_argument("--stage", choices=("info", "plot", "all"), default="all")
        if name in ("submit", "pipeline"):
            p.add_argument("--time", help="Slurm wall time")
        if name == "submit":
            p.add_argument("--test-only", action="store_true")
            p.add_argument("--resume", action="store_true", help="resubmit after previous job has left the queue")
    args = parser.parse_args()
    args.root = args.root.resolve()
    args.repo = (args.repo or args.root / "tools" / "dSph_workbench").resolve()
    try:
        return args.func(args)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

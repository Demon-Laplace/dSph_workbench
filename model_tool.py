#!/usr/bin/env python3
"""Create, edit, validate, and build models from IC INI files.

Examples
--------
Create a model from an existing model and change parameters::

    python3 model_tool.py create 5010 --base 4001 \
        --set Density_stellar.axisRatioZ=0.7 \
        --set Density_stellar.mass=0.008

Edit an existing model (linked component ratios are updated automatically)::

    python3 model_tool.py set 5010 Density_stellar.axisRatioZ=0.65

Check or repair existing files::

    python3 model_tool.py check "4002-4009,5001-5005"
    python3 model_tool.py sync "4002-4009,5001-5005"

Build model directories, replacing the execution role of makeini.py::

    python3 model_tool.py build "4002-4009,5001-5005" --processes 4

For compatibility, omitting ``build`` also works::

    python3 model_tool.py 5001

Before every write or build, ``axisRatioZ`` in each ``Component_*`` section
that has such a constraint is copied from the referenced ``Density_*``
section. Section and option names are matched case-insensitively. Before builds,
changed ``Fornax2047+`` parameter files are synchronized to ``origin/main``.
Successful builds automatically update and push ``IC_MODELS.csv``. Configure
``ic.root`` and ``ic_python`` in workflow.json; outputs use unique per-model ic/ runs.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime, timezone
import uuid


PROJECT_ROOT = Path(__file__).resolve().parent
INI_DIR = PROJECT_ROOT / "IniFile"
PROJECT_NAME = "Fornax"
DEFAULT_THREADS = 6
COMMANDS = {"create", "set", "check", "sync", "build"}
SYNC_MIN_MODEL = 2047


def sync_generated_models(paths):
    """Commit and push successfully generated Fornax2047+ IC files."""
    prefix = "IC_{}".format(PROJECT_NAME)
    eligible = []
    for path in paths:
        if path.name.startswith(prefix) and path.suffix == ".ini":
            number = path.stem[len(prefix):]
            if number.isdigit() and int(number) >= SYNC_MIN_MODEL:
                eligible.append(path)
    if not eligible:
        print("no eligible IC files to push")
        return
    git = ["git", "-C", str(PROJECT_ROOT)]
    try:
        subprocess.run(git + ["fetch", "origin", "main"], check=True)
        subprocess.run(git + ["merge", "--ff-only", "origin/main"], check=True)
        relative_paths = [str(path.relative_to(PROJECT_ROOT)) for path in eligible]
        subprocess.run(git + ["add", "-f"] + relative_paths, check=True)
        staged = subprocess.run(
            git + ["diff", "--cached", "--quiet", "--"] + relative_paths,
            check=False,
        )
        if staged.returncode == 0:
            print("no IC changes to push")
            return
        if staged.returncode != 1:
            raise subprocess.CalledProcessError(staged.returncode, staged.args)
        names = ", ".join(path.stem.removeprefix("IC_") for path in eligible)
        subprocess.run(
            git
            + ["commit", "--only", "-m", "sync {} IC models".format(names), "--"]
            + relative_paths,
            check=True,
        )
        subprocess.run(git + ["push", "origin", "main"], check=True)
        print("pushed {} IC file(s) to main".format(len(eligible)))
    except subprocess.CalledProcessError as error:
        raise IniError(
            "automatic Git sync failed with exit code {}".format(error.returncode)
        )


SECTION_RE = re.compile(r"^\s*\[([^]]+)]")
OPTION_RE = re.compile(
    r"^(\s*)([^#;\s][^=:]*?)(\s*)([=:])(\s*)(.*?)(\r?\n?)$"
)


class IniError(ValueError):
    """Raised for an invalid or unsafe INI edit."""


class IniDocument:
    """Small line-preserving INI editor for the project's simple INI files."""

    def __init__(self, text, source="<memory>"):
        self.lines = text.splitlines(keepends=True)
        self.source = str(source)
        self._scan()

    @classmethod
    def read(cls, path):
        path = Path(path)
        return cls(path.read_text(encoding="utf-8"), path)

    def _scan(self):
        self.sections = {}
        self.options = {}
        current = None
        for index, line in enumerate(self.lines):
            section_match = SECTION_RE.match(line)
            if section_match:
                name = section_match.group(1).strip()
                normalized = name.casefold()
                if normalized in self.sections:
                    raise IniError(
                        "{}: duplicate section [{}]".format(self.source, name)
                    )
                self.sections[normalized] = name
                self.options[normalized] = {}
                current = normalized
                continue
            if current is None:
                continue
            option_match = OPTION_RE.match(line)
            if not option_match:
                continue
            key = option_match.group(2).strip()
            normalized_key = key.casefold()
            if normalized_key in self.options[current]:
                raise IniError(
                    "{}: duplicate option {} in [{}]".format(
                        self.source, key, self.sections[current]
                    )
                )
            self.options[current][normalized_key] = (index, key, option_match)

    def section_names(self):
        return list(self.sections.values())

    def has_option(self, section, option):
        normalized = section.casefold()
        return (
            normalized in self.options
            and option.casefold() in self.options[normalized]
        )

    def get(self, section, option):
        normalized = section.casefold()
        try:
            _, _, match = self.options[normalized][option.casefold()]
        except KeyError:
            raise IniError(
                "{}: missing {}.{}".format(self.source, section, option)
            )
        return match.group(6).strip()

    def set_existing(self, section, option, value):
        normalized = section.casefold()
        normalized_option = option.casefold()
        try:
            index, _, match = self.options[normalized][normalized_option]
        except KeyError:
            raise IniError(
                "{}: refusing to create unknown setting {}.{}; check the spelling"
                .format(self.source, section, option)
            )
        self.lines[index] = "{}{}{}{}{}{}{}".format(
            match.group(1),
            match.group(2),
            match.group(3),
            match.group(4),
            match.group(5),
            value,
            match.group(7),
        )
        self._scan()

    def text(self):
        return "".join(self.lines)


def values_equal(first, second):
    """Treat numerically equivalent spellings (for example 0.6 and 0.60) alike."""
    try:
        return float(first) == float(second)
    except ValueError:
        return first == second


def ratio_links(document):
    """Return linked component/density ratios and structural errors."""
    links = []
    errors = []
    for component in document.section_names():
        if not component.casefold().startswith("component_"):
            continue
        if not document.has_option(component, "axisRatioZ"):
            continue
        if not document.has_option(component, "density"):
            errors.append("[{}] has axisRatioZ but no density reference".format(component))
            continue
        density = document.get(component, "density")
        if density.casefold() not in document.sections:
            errors.append(
                "[{}] references missing section [{}]".format(component, density)
            )
            continue
        density_name = document.sections[density.casefold()]
        if not document.has_option(density_name, "axisRatioZ"):
            errors.append(
                "[{}] has axisRatioZ but [{}] does not".format(
                    component, density_name
                )
            )
            continue
        links.append(
            (
                component,
                density_name,
                document.get(component, "axisRatioZ"),
                document.get(density_name, "axisRatioZ"),
            )
        )
    return links, errors


def synchronize_ratios(document):
    """Copy linked density ratios into component constraints."""
    links, errors = ratio_links(document)
    if errors:
        raise IniError("{}: {}".format(document.source, "; ".join(errors)))
    changes = []
    for component, density, component_value, density_value in links:
        if not values_equal(component_value, density_value):
            document.set_existing(component, "axisRatioZ", density_value)
            changes.append(
                "[{}].axisRatioZ: {} -> {} (from [{}])".format(
                    component, component_value, density_value, density
                )
            )
    return changes


def ratio_problems(document):
    links, errors = ratio_links(document)
    problems = list(errors)
    for component, density, component_value, density_value in links:
        if not values_equal(component_value, density_value):
            problems.append(
                "[{}].axisRatioZ={} but [{}].axisRatioZ={}".format(
                    component, component_value, density, density_value
                )
            )
    return problems


def log_tail(path, max_lines=40, max_characters=8000):
    """Read a bounded diagnostic tail without hiding the full log location."""
    try:
        lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as error:
        return "unable to read log: {}".format(error)
    text = "\n".join(lines[-max_lines:])
    if len(text) > max_characters:
        text = "...\n" + text[-max_characters:]
    return text or "log is empty"


def write_atomic(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = path.stat().st_mode if path.exists() else None
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=str(path.parent), delete=False
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    if mode is not None:
        os.chmod(str(temporary), mode)
    os.replace(str(temporary), str(path))


def parse_overrides(items):
    overrides = []
    for item in items:
        if "=" not in item or "." not in item.split("=", 1)[0]:
            raise IniError(
                "invalid --set {!r}; expected SECTION.OPTION=VALUE".format(item)
            )
        target, value = item.split("=", 1)
        section, option = target.rsplit(".", 1)
        if not section or not option or not value:
            raise IniError(
                "invalid --set {!r}; expected SECTION.OPTION=VALUE".format(item)
            )
        overrides.append((section, option, value))
    return overrides


def model_ini_path(value):
    candidate = Path(value)
    if candidate.suffix.casefold() == ".ini" or candidate.parent != Path("."):
        return candidate.resolve()
    name = value
    if name.startswith("IC_"):
        filename = name + ".ini"
    elif name.startswith(PROJECT_NAME):
        filename = "IC_" + name + ".ini"
    elif name.isdigit():
        filename = "IC_{}{}.ini".format(PROJECT_NAME, name)
    else:
        raise IniError("cannot interpret model or INI path {!r}".format(value))
    return INI_DIR / filename


def expand_model_arguments(arguments):
    if len(arguments) != 1 or not re.fullmatch(r"[\d,\-\s]+", arguments[0]):
        return arguments
    numbers = set()
    for item in arguments[0].replace(" ", "").split(","):
        if not item:
            raise IniError("empty item in model selection")
        if "-" in item:
            parts = item.split("-")
            if len(parts) != 2 or not all(part.isdigit() for part in parts):
                raise IniError("invalid model range {!r}".format(item))
            start, end = (int(part) for part in parts)
            if start > end:
                raise IniError("descending model range {!r}".format(item))
            numbers.update(range(start, end + 1))
        elif item.isdigit():
            numbers.add(int(item))
        else:
            raise IniError("invalid model number {!r}".format(item))
    return [str(number) for number in sorted(numbers)]


def existing_ini_paths(arguments):
    paths = [model_ini_path(item) for item in expand_model_arguments(arguments)]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise IniError("missing INI file(s): {}".format(", ".join(missing)))
    return paths


def available_cpu_count():
    """Return CPUs available to this process, respecting scheduler affinity."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def command_create(args):
    source = model_ini_path(args.base)
    target = model_ini_path(args.model)
    if not source.is_file():
        raise IniError("base INI does not exist: {}".format(source))
    if target.exists() and not args.force:
        raise IniError("target already exists (use --force to replace): {}".format(target))
    document = IniDocument.read(source)
    document.source = str(target)
    for section, option, value in parse_overrides(args.settings):
        document.set_existing(section, option, value)
    changes = synchronize_ratios(document)
    write_atomic(target, document.text())
    print("created {}".format(target))
    for change in changes:
        print("  synced {}".format(change))
    return 0


def command_set(args):
    path = model_ini_path(args.model)
    if not path.is_file():
        raise IniError("INI file does not exist: {}".format(path))
    document = IniDocument.read(path)
    for section, option, value in parse_overrides(args.settings):
        document.set_existing(section, option, value)
    changes = synchronize_ratios(document)
    write_atomic(path, document.text())
    print("updated {}".format(path))
    for change in changes:
        print("  synced {}".format(change))
    return 0


def command_check(args):
    failed = False
    for path in existing_ini_paths(args.models):
        problems = ratio_problems(IniDocument.read(path))
        if problems:
            failed = True
            print("ERROR {}".format(path))
            for problem in problems:
                print("  {}".format(problem))
        else:
            print("OK {}".format(path))
    return 1 if failed else 0


def command_sync(args):
    for path in existing_ini_paths(args.models):
        document = IniDocument.read(path)
        changes = synchronize_ratios(document)
        if changes:
            write_atomic(path, document.text())
            print("updated {}".format(path))
            for change in changes:
                print("  {}".format(change))
        else:
            print("unchanged {}".format(path))
    return 0


def build_one(task):
    model_name, ini_path, threads, seed = task
    name = model_name + '-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
    command = [sys.executable, str(PROJECT_ROOT / 'workflow/workflow.py'),
               '--root', str(PROJECT_ROOT), '--repo', str(PROJECT_ROOT),
               'build-ic', '--ini', str(ini_path), '--name', name,
               '--threads', str(threads), '--seed', str(seed)]
    result = subprocess.run(command)
    return model_name, result.returncode == 0, f'Run {name}; exit status {result.returncode}'


def command_build(args):
    if args.overwrite:
        raise IniError("--overwrite is retired; each build now creates a separate per-model ic/ run")
    paths = existing_ini_paths(args.models)
    models = []
    for path in paths:
        document = IniDocument.read(path)
        changes = synchronize_ratios(document)
        if changes:
            write_atomic(path, document.text())
            print("preflight synced {}".format(path))
            for change in changes:
                print("  {}".format(change))
        stem = path.stem
        model_name = stem[3:] if stem.startswith("IC_") else stem
        models.append((model_name, path))

    if args.processes is not None and args.processes < 1:
        raise IniError("--processes must be at least 1")
    if args.threads is not None and args.threads < 1:
        raise IniError("--threads must be at least 1")
    if not models:
        print("no models selected for building")
        return 0

    available = available_cpu_count()
    requested_processes = args.processes or min(
        len(models), max(1, available // DEFAULT_THREADS)
    )
    # Processes parallelize models, so more workers than models only add overhead.
    processes = min(requested_processes, len(models))
    if requested_processes != processes:
        print(
            "using {} process(es), because only {} model(s) were requested".format(
                processes, len(models)
            )
        )

    threads = args.threads or min(
        DEFAULT_THREADS, max(1, available // processes)
    )
    if processes * threads > available:
        print(
            "WARNING: {} processes x {} threads exceeds {} available CPUs".format(
                processes, threads, available
            )
        )
    # Preserve the established parameter synchronization policy before builds.
    sync_generated_models(paths)
    tasks = [(name, path, threads, args.seed) for name, path in models]
    print(
        "building {} model(s) with {} process(es) and {} AGAMA thread(s) per model"
        .format(len(tasks), processes, threads)
    )

    results = []
    if processes == 1:
        for task in tasks:
            results.append(build_one(task))
    else:
        with ProcessPoolExecutor(max_workers=processes) as executor:
            futures = [executor.submit(build_one, task) for task in tasks]
            for future in as_completed(futures):
                results.append(future.result())

    failures = []
    for model_name, success, error in sorted(results):
        if success:
            print("OK {}".format(model_name))
        else:
            failures.append((model_name, error))
            print("ERROR {}\n{}".format(model_name, error))
    return 1 if failures else 0


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    create = subparsers.add_parser("create", help="copy a base INI and apply settings")
    create.add_argument("model", help="new model number, name, or INI path")
    create.add_argument("--base", required=True, help="base model number, name, or INI path")
    create.add_argument(
        "--set", dest="settings", action="append", default=[],
        metavar="SECTION.OPTION=VALUE", help="setting to override; may be repeated"
    )
    create.add_argument("--force", action="store_true", help="replace an existing target INI")
    create.set_defaults(handler=command_create)

    edit = subparsers.add_parser("set", help="edit settings and synchronize linked ratios")
    edit.add_argument("model", help="model number, name, or INI path")
    edit.add_argument("settings", nargs="+", metavar="SECTION.OPTION=VALUE")
    edit.set_defaults(handler=command_set)

    check = subparsers.add_parser("check", help="report linked ratio inconsistencies")
    check.add_argument("models", nargs="+", help="model selection or INI paths")
    check.set_defaults(handler=command_check)

    sync = subparsers.add_parser("sync", help="repair linked ratio inconsistencies")
    sync.add_argument("models", nargs="+", help="model selection or INI paths")
    sync.set_defaults(handler=command_sync)

    build = subparsers.add_parser("build", help="generate model directories from INI files")
    build.add_argument("models", nargs="+", help="model selection or INI paths")
    build.add_argument(
        "--processes", type=int,
        help="parallel model count (does not accelerate one model)"
    )
    build.add_argument(
        "--threads", type=int,
        help="OpenMP threads per AGAMA model (default: 6, capped by available CPUs)"
    )
    build.add_argument(
        "--overwrite", action="store_true",
        help="retired; existing outputs are always preserved"
    )
    build.add_argument("--seed", type=int, default=42)
    build.set_defaults(handler=command_build)
    return parser


def main(argv=None):
    parser = make_parser()
    arguments = list(sys.argv[1:] if argv is None else argv)
    if (
        arguments
        and arguments[0] not in COMMANDS
        and arguments[0] not in ("-h", "--help")
    ):
        arguments.insert(0, "build")
    args = parser.parse_args(arguments)
    try:
        return args.handler(args)
    except IniError as error:
        parser.error(str(error))


if __name__ == "__main__":
    sys.exit(main())

"""Small run-version ledger; never opens particle inputs or computes hashes."""
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import tempfile
from datetime import datetime, timezone


CONTRACT = 'isolation_run_versions_v1'
LEDGER = 'MODEL_VERSIONS.json'
TEXT = 'MODEL_VERSIONS.txt'


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_text(path, text):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile('w', encoding='utf-8', dir=path.parent,
                                         prefix='.' + path.name + '.', delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(text)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def source(request):
    """A declared build/version, not a claim about the bytes in its files."""
    build = request.get('build_id')
    if build:
        result = {'build_id': str(build)}
        if request.get('source_version'):
            result['declared_version'] = str(request['source_version'])
        return result
    if request.get('source_version'):
        return {'declared_version': str(request['source_version'])}
    directory = request.get('source_directory')
    return {'directory': str(directory)} if directory else None


def same_source(saved, requested):
    a, b = source(saved), source(requested)
    if saved.get('source_version') != requested.get('source_version'):
        return False
    if a and b and 'build_id' in a and 'build_id' in b:
        return a == b
    if bool(saved.get('build_id')) != bool(requested.get('build_id')):
        return False
    if saved.get('source_version'):
        return a == b
    a, b = saved.get('source_directory'), requested.get('source_directory')
    return bool(a and b and a == b)


def parameter_values(text):
    values = {}
    for line in text.splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2:
            values[fields[0]] = fields[1]
    # These are generated for each run, not scientific configuration changes.
    for key in ('InitCondFile', 'OutputDir', 'ResubmitOn'):
        values.pop(key, None)
    return values


def parameter_signature(profile, settings, template_text):
    values = parameter_values(template_text)
    values.update({key: str(value) for key, value in settings.get('parameter_overrides', {}).items()
                   if key not in ('InitCondFile', 'OutputDir', 'ResubmitOn')})
    configuration = {key: value for key, value in settings.items()
                     if key != 'parameter_overrides'}
    return {'profile': profile, 'parameters': values, 'configuration': configuration}


def empty(model):
    return {'contract': CONTRACT, 'model': model, 'ics': {}, 'parameters': {}, 'runs': {}}


@contextmanager
def locked(directory, model):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name in (LEDGER, TEXT, '.model-versions.lock'):
        path = directory / name
        if path.is_symlink():
            raise ValueError('Version bookkeeping must stay in this model directory: ' + str(path))
    with (directory / '.model-versions.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        path = directory / LEDGER
        data = json.loads(path.read_text()) if path.exists() else empty(model)
        if data.get('contract') != CONTRACT or data.get('model') != model:
            raise ValueError('Version ledger belongs to another model or schema')
        yield data


def save(directory, data):
    """JSON is authoritative. TXT is a replaceable view of the same records."""
    atomic_text(directory / LEDGER, json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    lines = ['Model: ' + data['model'],
             'Declared IC/configuration versions; no particle or content verification.',
             'State complete means solver process completion, not equilibrium or scientific validation.',
             'JSON source: ' + LEDGER, '']
    for name, entry in data['runs'].items():
        request = entry['request']
        lines.extend([
            '[' + name + ']',
            'requested_name: ' + request.get('requested_name', name),
            'ic_version: ' + str(request.get('ic_version') or 'legacy-unversioned'),
            'parameter_version: ' + str(request.get('parameter_version') or 'legacy-unversioned'),
            'build_id: ' + str(request.get('build_id') or request.get('source_version') or 'unknown'),
            'source_directory: ' + str(request.get('source_directory') or 'unknown'),
            'profile: ' + str(request.get('profile') or 'unknown'),
            'state: ' + entry.get('state', 'reserved'),
            'job_id: ' + str(entry.get('job_id') or '-'),
            'directory: ' + entry['directory'], '',
        ])
    try:
        atomic_text(directory / TEXT, '\n'.join(lines))
    except OSError as error:
        # A stale/missing human-readable view must not invalidate a job binding.
        import warnings
        warnings.warn('Version JSON saved; TXT view could not be refreshed: ' + str(error))


def number(mapping, prefix, value, key=None):
    for version, record in mapping.items():
        if (record[key] if key else record) == value:
            return version
    version = prefix + str(max([int(name[len(prefix):]) for name in mapping] or [0]) + 1).zfill(3)
    mapping[version] = value if key is None else {key: value}
    return version


def resolve(directory, request, settings, template_text, locate):
    """Reserve one stable run name per requested name and IC/parameter version."""
    directory = Path(directory)
    request = dict(request)
    model = request['model']
    requested_name = request.get('requested_name', request['name'])
    original = locate(requested_name)
    signature = parameter_signature(request['profile'], settings, template_text)
    declared_source = source(request)
    if declared_source is None:
        # Untagged imported inputs are tied to their explicitly supplied location.
        declared_source = {'directory': request.get('source_directory', 'legacy-unspecified')}
    with locked(directory, model) as data:
        changed = False
        # Register, but never rewrite or relabel, an existing historical run.
        if requested_name not in data['runs'] and original.exists():
            old = json.loads((original / 'request.json').read_text())
            if old.get('model') != model or old.get('name') != requested_name:
                raise ValueError('Existing run directory belongs to another request')
            old = dict(old, requested_name=requested_name, version_contract='legacy')
            state_path = original / 'status.json'
            status = json.loads(state_path.read_text()) if state_path.exists() else {}
            data['runs'][requested_name] = {
                'request': old, 'directory': str(original),
                'state': status.get('state', 'legacy'), 'job_id': status.get('job_id')}
            profile_path = original / 'profile.json'
            parameter_path = original / 'GZWJL.PARAM'
            if profile_path.is_file() and parameter_path.is_file():
                data['runs'][requested_name]['legacy_parameters'] = parameter_signature(
                    old.get('profile'), json.loads(profile_path.read_text()), parameter_path.read_text())
            changed = True
        # Explicit actual run names identify that run, not an unrelated alias.
        direct = data['runs'].get(request['name'])
        if direct:
            requested_name = direct['request'].get('requested_name', request['name'])
        for entry in data['runs'].values():
            old = entry['request']
            if old.get('requested_name', old['name']) != requested_name:
                continue
            if old.get('version_contract') == 'legacy':
                compatible = (old.get('profile') == request['profile'] and same_source(old, request)
                              and ('legacy_parameters' not in entry
                                   or entry['legacy_parameters'] == signature))
                if compatible:
                    if changed:
                        save(directory, data)
                    return dict(entry, reused=True)
                continue
            if (data['ics'][old['ic_version']] == declared_source
                    and data['parameters'][old['parameter_version']]['signature'] == signature):
                if changed:
                    save(directory, data)
                return dict(entry, reused=True)
        old_parameters = len(data['parameters'])
        ic = number(data['ics'], 'IC', declared_source)
        parameter = number(data['parameters'], 'P', signature, key='signature')
        if len(data['parameters']) > old_parameters:
            data['parameters'][parameter].update(settings=settings, template_text=template_text)
        name = requested_name
        if name in data['runs'] or locate(name).exists():
            stem = requested_name + '-' + ic + '-' + parameter
            suffix = 1
            name = stem + '-' + str(suffix).zfill(2)
            while name in data['runs'] or locate(name).exists():
                suffix += 1
                name = stem + '-' + str(suffix).zfill(2)
        request.update(name=name, requested_name=requested_name, version_contract=CONTRACT,
                       ic_version=ic, parameter_version=parameter)
        entry = {'request': request, 'directory': str(locate(name)),
                 'state': 'reserved', 'created_utc': now()}
        data['runs'][name] = entry
        save(directory, data)
        return dict(entry, reused=False)


def binding(directory, request):
    """Use the configuration snapshot already selected by resolve, not live settings."""
    with locked(directory, request['model']) as data:
        entry = data['runs'].get(request['name'])
        keys = ('model', 'name', 'profile', 'ic_version', 'parameter_version')
        if (not entry or any(entry['request'].get(key) != request.get(key) for key in keys)
                or not same_source(entry['request'], request)):
            raise ValueError('Requested run/version binding is not registered')
        parameters = data['parameters'][request['parameter_version']]
        return dict(entry), parameters['settings'], parameters['template_text']


def lookup(directory, model, name):
    directory = Path(directory)
    if not (directory / LEDGER).is_file():
        return None
    with locked(directory, model) as data:
        return data['runs'].get(name)


def update(directory, model, name, status):
    directory = Path(directory)
    if not (directory / LEDGER).is_file():
        return
    with locked(directory, model) as data:
        entry = data['runs'].get(name)
        if entry is None:
            return
        fields = {key: status[key] for key in ('state', 'job_id') if key in status}
        if all(entry.get(key) == value for key, value in fields.items()):
            return
        entry.update(fields, updated_utc=now())
        save(directory, data)


def continuation(directory, request, settings, text, previous):
    """Keep one physical run, preserving the earlier binding as attempt history."""
    directory = Path(directory)
    with locked(directory, request['model']) as data:
        entry = data['runs'].setdefault(request['name'], {
            'request': dict(request, version_contract=request.get('version_contract', 'legacy')),
            'directory': str(directory / 'runs' / request['name'])})
        signature = parameter_signature(request.get('profile'), settings, text)
        parameter = number(data['parameters'], 'P', signature, key='signature')
        data['parameters'][parameter].update(settings=settings, template_text=text)
        entry.setdefault('continuations', []).append(
            {'request': dict(entry['request']), 'status': previous, 'recorded_utc': now()})
        entry['request'] = dict(entry['request'], parameter_version=parameter)
        entry.update(state='submitting', updated_utc=now())
        save(directory, data)
        return parameter

"""Bind production TAN media to its display and physical input identities."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

from canonical_display import TAN_DISPLAY_CONTRACT_ID, DISPLAY_RELATIVE_DIR


RENDER_REVISION = 5
DISPLAY_IDENTITY_FIELDS = (
    'canonical_observer_id', 'member_selection_id', 'display_contract_id',
    'relative_directory', 'render_revision',
)
_ordinary_display_identities = {}


def ordinary_media_reuse():
    """Explicit job mode: reuse the same scientific definition without media QA."""
    return os.environ.get('DSPH_MEDIA_REUSE') == 'scientific'


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write('\n')
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def display_identity(summary):
    identities = []
    with Path(summary).open(encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('# canonical_observer_id='):
                identities.append(dict(item.split('=', 1) for item in line[2:].split() if '=' in item))
            elif line.strip() and not line.startswith('#'):
                break
    if len(identities) != 1:
        raise ValueError('Display requires one canonical summary identity')
    identity = identities[0]
    required = ('canonical_observer_id', 'member_selection_id')
    if any(not identity.get(key) for key in required):
        raise ValueError('Incomplete canonical summary identity')
    if ordinary_media_reuse():
        from canonical_observer import require_canonical_csv
        require_canonical_csv(summary)
        identity.update(display_contract_id=TAN_DISPLAY_CONTRACT_ID,
                        relative_directory=DISPLAY_RELATIVE_DIR,
                        render_revision=RENDER_REVISION)
        return {key: identity[key] for key in DISPLAY_IDENTITY_FIELDS}
    source = Path(__file__).resolve().parent
    identity.update(display_contract_id=TAN_DISPLAY_CONTRACT_ID,
                    display_module_sha256=file_sha256(source / 'canonical_display.py'),
                    relative_directory=DISPLAY_RELATIVE_DIR,
                    render_source_sha256={path.name: file_sha256(path)
                                          for path in sorted(source.glob('*.py'))})
    return identity


def frame_is_current(path, identity):
    path = Path(path)
    sidecar = Path(str(path) + '.display.json')
    if ordinary_media_reuse():
        if not path.is_file():
            return False
        manifest = path.parent.parent / 'display_manifest.json'
        key = str(manifest.resolve())
        if key not in _ordinary_display_identities:
            if not manifest.is_file():
                return False
            _ordinary_display_identities[key] = json.loads(manifest.read_text(encoding='utf-8'))
        record = _ordinary_display_identities[key]
        return all(record.get(key) == identity.get(key)
                   for key in DISPLAY_IDENTITY_FIELDS)
    if not path.is_file() or not sidecar.is_file():
        return False
    try:
        record = json.loads(sidecar.read_text(encoding='utf-8'))
        return all(record.get(key) == value for key, value in identity.items()) and record.get('png_sha256') == file_sha256(path)
    except (ValueError, OSError):
        return False


def record_frame(path, identity, numsp):
    path = Path(path)
    if ordinary_media_reuse():
        return
    write_json_atomic(str(path) + '.display.json', dict(identity, numsp=int(numsp),
                      png_sha256=file_sha256(path), bytes=path.stat().st_size))


def record_display_manifest(display_root, identity, snapshots, video_paths=()):
    root = Path(display_root)
    if ordinary_media_reuse():
        videos = []
        for path in video_paths:
            path = Path(path)
            if path.parent.resolve() != root.resolve():
                raise ValueError('Video is outside this display output directory')
            videos.append(dict(name=path.name))
        write_json_atomic(root / 'display_manifest.json', dict(
            identity, frames=[dict(numsp=int(n)) for n in snapshots], videos=videos))
        _ordinary_display_identities[str((root / 'display_manifest.json').resolve())] = identity
        return
    frames = []
    for number in snapshots:
        path = root / 'output_img' / ('check_nature_%03d.png' % number)
        if not frame_is_current(path, identity):
            raise ValueError('Missing or different-display frame: ' + str(path))
        frames.append(dict(numsp=int(number), sha256=file_sha256(path), bytes=path.stat().st_size))
    videos = []
    for path in video_paths:
        path = Path(path)
        if path.parent.resolve() != root.resolve() or not path.is_file():
            raise ValueError('Video is outside this display output directory')
        videos.append(dict(name=path.name, sha256=file_sha256(path), bytes=path.stat().st_size))
    write_json_atomic(root / 'display_manifest.json', dict(identity, frames=frames, videos=videos))


def prepare_display_directory(display_root, identity):
    """Start/resume this explicitly scoped scientific output, without frame QA."""
    if not ordinary_media_reuse():
        return
    from canonical_observer import MEMBER_SELECTION_ID
    root = Path(display_root)
    analysis_root = root.parents[1]
    if analysis_root.name != MEMBER_SELECTION_ID:
        raise ValueError('Ordinary media reuse requires a dedicated ' + MEMBER_SELECTION_ID + ' analysis root')
    manifest = root / 'display_manifest.json'
    if manifest.is_file():
        record = json.loads(manifest.read_text(encoding='utf-8'))
        if any(record.get(key) != identity.get(key) for key in DISPLAY_IDENTITY_FIELDS):
            raise ValueError('Different scientific display definition; use a fresh analysis root: ' + str(root))
        _ordinary_display_identities[str(manifest.resolve())] = record
    else:
        # The member-named root isolates this batch from all historical frames.
        record_display_manifest(root, identity, [])

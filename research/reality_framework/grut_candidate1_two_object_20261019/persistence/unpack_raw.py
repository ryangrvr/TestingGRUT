"""Validate and unpack only the manifest's raw diagnostic files, preserving changes."""
from pathlib import Path, PurePosixPath
import hashlib
import io
import json
import tarfile

ROOT = Path(__file__).resolve().parent


def allowed(name):
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or str(path) != name:
        return False
    parts = list(path.parts)
    if parts and parts[0] in ('formation', 'record'):
        parts.pop(0)
    return (len(parts) == 1 and parts[0] == 'execution_status.json') or (
        len(parts) == 2 and ((parts[0] == 'results' and parts[1].endswith('.json')) or
                             (parts[0] == 'logs' and parts[1].endswith('.log'))))


def main():
    manifest = json.loads((ROOT / 'RAW_OUTPUTS_MANIFEST.json').read_text())
    if manifest['archive'] != 'raw_outputs.tar.gz':
        raise SystemExit('Unexpected archive path')
    archive = (ROOT / manifest['archive']).read_bytes()
    if len(archive) != manifest['bytes'] or hashlib.sha256(archive).hexdigest() != manifest['sha256']:
        raise SystemExit('Archive identity mismatch')
    expected = {x['path']: x for x in manifest['entries']}
    if len(expected) != len(manifest['entries']):
        raise SystemExit('Duplicate manifest path')
    payloads = {}
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r:gz') as tar:
        for member in tar:
            name = member.name
            if not member.isfile() or not allowed(name) or name not in expected or name in payloads:
                raise SystemExit('Unapproved archive entry: ' + name)
            data = tar.extractfile(member).read()
            entry = expected[name]
            if len(data) != entry['bytes'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
                raise SystemExit('Entry identity mismatch: ' + name)
            dest = ROOT / name
            if not dest.resolve().is_relative_to(ROOT.resolve()):
                raise SystemExit('Destination escapes package: ' + name)
            if dest.is_symlink() or (dest.exists() and (not dest.is_file() or dest.read_bytes() != data)):
                raise SystemExit('Existing file differs; preserve or move it before unpacking: ' + name)
            payloads[name] = data
    if set(payloads) != set(expected):
        raise SystemExit('Archive entry set differs from manifest')
    # Validation completes before the first write. Repeated identical unpacking is harmless.
    for name, data in payloads.items():
        dest = ROOT / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            dest.write_bytes(data)
    print(json.dumps(dict(verified_files=len(payloads), archive_sha256=manifest['sha256'])))


if __name__ == '__main__':
    main()

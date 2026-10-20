"""Package this execution's raw diagnostics with deterministic archive metadata."""
from pathlib import Path
import gzip
import hashlib
import io
import json
import tarfile

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / 'raw_outputs.tar.gz'


def main():
    paths = []
    for study in ['', 'formation', 'record']:
        folder = ROOT / study
        paths.extend(sorted((folder / 'results').glob('*.json')))
        paths.extend(sorted((folder / 'logs').glob('*.log')))
        paths.append(folder / 'execution_status.json')
    if sum(p.parent.name == 'results' for p in paths) != 31:
        raise SystemExit('Expected exactly 31 completed result files')
    entries = []
    with ARCHIVE.open('wb') as out:
        with gzip.GzipFile(filename='', mode='wb', fileobj=out, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode='w|', format=tarfile.USTAR_FORMAT) as tar:
                for path in sorted(paths):
                    data = path.read_bytes()
                    name = path.relative_to(ROOT).as_posix()
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = 0o644
                    info.mtime = 0
                    tar.addfile(info, io.BytesIO(data))
                    entries.append(dict(path=name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    result = dict(archive=ARCHIVE.name, bytes=ARCHIVE.stat().st_size,
                  sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(), entries=entries,
                  note='Exact local execution outputs. Timestamps inside logs and wall_seconds are execution metadata; archive metadata is deterministic.')
    (ROOT / 'RAW_OUTPUTS_MANIFEST.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(files=len(entries), bytes=result['bytes'], sha256=result['sha256'])))


if __name__ == '__main__':
    main()

"""Archive the ten exact full-field diagnostic histories and execution status."""
from pathlib import Path
import gzip
import hashlib
import io
import json
import tarfile

ROOT=Path(__file__).resolve().parent

if __name__=='__main__':
    protocol=json.loads((ROOT/'protocol.json').read_text())
    paths=[ROOT/'results'/(spec['id']+'.json') for spec in protocol['runs']]+[ROOT/'execution_status.json']
    archive=ROOT/'raw_outputs.tar.gz'
    entries=[]
    with archive.open('wb') as output:
        with gzip.GzipFile(filename='',mode='wb',fileobj=output,mtime=0) as gz:
            with tarfile.open(fileobj=gz,mode='w|',format=tarfile.USTAR_FORMAT) as tar:
                for path in sorted(paths):
                    data=path.read_bytes();name=path.relative_to(ROOT).as_posix()
                    info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644;info.mtime=0
                    tar.addfile(info,io.BytesIO(data))
                    entries.append(dict(path=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    result=dict(archive=archive.name,bytes=archive.stat().st_size,
                sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),entries=entries)
    (ROOT/'RAW_OUTPUTS_MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(files=len(entries),bytes=result['bytes'],sha256=result['sha256'])))

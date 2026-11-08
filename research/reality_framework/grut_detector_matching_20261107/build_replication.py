"""Build a small standalone source/results archive from a clean committed tree."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    if subprocess.check_output(["git","status","--porcelain"],cwd=REPO,text=True):
        raise RuntimeError("Package must be built from a clean committed tree")
    head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=REPO,text=True).strip()
    prefixes=[str(HERE.relative_to(REPO)),"research/reality_framework/grut_qed_detector_histories_20261107"]
    files=subprocess.check_output(["git","ls-files","-z","--",*prefixes],cwd=REPO).decode().split("\0")
    files=[name for name in files if name and not name.endswith(".zip")]
    manifest={"archive_source_commit":head,"physics_freezes":{"M0-01":"2b13287ae398d45bed7b4dcdcbb083747aaa4a00","M0-02":"c5e107a9920c57b13482d63f290fa9f7bc1fc7b5","M1-01":"479eaed2f906e323dac491f8f308c0c05718d63b","D1-01":"90745ef0d7aaf4c6d8fe3a11edc7ee0d12f7d339"},
              "sha256":{name:hashlib.sha256((REPO/name).read_bytes()).hexdigest() for name in sorted(files)}}
    archive_readme='''# GRUT interacting-history replication packet

This archive contains both independent-history source/results and the original QED packet it imports. It is not a publication release. All matching/reduction misses are retained. Zenodo remains paused.

Start with `research/reality_framework/grut_detector_matching_20261107/README.md`. The manifest records the original source commit and physics freezes. Extraction does not recreate historical Git metadata. The documented local bootstrap creates a new commit solely to allow the clean-tree execution guards; source hashes retain the historical comparison.

From the extracted archive root, before installing or executing anything:

```sh
python - <<'PY'
import hashlib,json
from pathlib import Path
m=json.loads(Path('PACKET_MANIFEST.json').read_text())
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in m['sha256'].items())
print('Manifest matches',len(m['sha256']),'files; source',m['archive_source_commit'])
PY
git init -q
git add .
git -c user.name='Local replica' -c user.email='replica@localhost' commit -q -m 'Import verified replication packet'
```

Then follow the repository instructions in the nested README. Put the environment and outputs outside the extracted tree. Only NumPy and SciPy are needed for the new reference. TeNPy is optional and used only for the original four-case DMRG replay. Installing these dependencies may require network access.
'''
    with zipfile.ZipFile(args.output,"x",compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for name in sorted(files):z.write(REPO/name,name)
        z.writestr("PACKET_MANIFEST.json",json.dumps(manifest,indent=2)+"\n")
        z.writestr("START_HERE.md",archive_readme)
    print(json.dumps({"output":str(args.output.resolve()),"files":len(files),"source_commit":head,
                      "bytes":args.output.stat().st_size,"sha256":hashlib.sha256(args.output.read_bytes()).hexdigest()}))


if __name__=="__main__":main()

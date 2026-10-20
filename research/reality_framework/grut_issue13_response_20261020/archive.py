#!/usr/bin/env python3
"""Deterministic archive of new diagnostics; exact endpoint deltas against parent.

The six rerun endpoint volumes differ from the parent by only 83 bytes in
total. Store those explicit byte replacements and both hashes, rather than
duplicate 10 MB of existing data. This reconstructs exact recorded bytes;
it does not round away the small rerun differences.
"""
from pathlib import Path
import gzip
import hashlib
import json
import sys
import tarfile

P=Path(__file__).resolve().parent
OLD=P.parent/'grut_3d_empirical_20261020/raw'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def pack():
    files=sorted((P/'runs').iterdir())
    members=[];patches={}
    for f in files:
        if f.suffix!='.f32':
            members.append(f);continue
        old=OLD/f.name
        a,b=old.read_bytes(),f.read_bytes()
        if len(a)!=len(b):raise RuntimeError('Endpoint length changed')
        patches[f.name]={'parent_file':f.name,'parent_sha256':sha(old),
                         'byte_replacements':[[i,y] for i,(x,y) in enumerate(zip(a,b)) if x!=y]}
    target=P/'RAW_DIAGNOSTICS.tar.gz'
    with target.open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0,compresslevel=9) as zipped:
            with tarfile.open(fileobj=zipped,mode='w') as tar:
                for f in members:
                    info=tarfile.TarInfo('runs/'+f.name);info.size=f.stat().st_size;info.mode=0o644;info.mtime=0
                    with f.open('rb') as data:tar.addfile(info,data)
    manifest={'archive':{'name':target.name,'sha256':sha(target),'bytes':target.stat().st_size},
              'files':{f.name:{'sha256':sha(f),'bytes':f.stat().st_size} for f in files},
              'endpoint_reconstruction':patches,
              'parent_packet_commit':'479f111117ea7162bd7888ff8f3792b2f6c45fdb',
              'scope':'30 exact new diagnostic files in tar plus six exact endpoint reconstructions from frozen parent volumes. No full evolving fields/momenta were saved.'}
    (P/'RAW_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Archived',len(members),'members;',sum(len(x['byte_replacements']) for x in patches.values()),'endpoint byte replacements;',target.stat().st_size,'compressed bytes')

def unpack():
    m=json.loads((P/'RAW_MANIFEST.json').read_text());f=P/m['archive']['name']
    if sha(f)!=m['archive']['sha256']:raise RuntimeError('Archive changed')
    with tarfile.open(f,'r:gz') as tar:
        expected={'runs/'+name for name in m['files'] if name not in m['endpoint_reconstruction']}
        if {x.name for x in tar.getmembers()}!=expected or not all(x.isfile() for x in tar.getmembers()):
            raise RuntimeError('Unexpected archive members')
        tar.extractall(P,filter='data')
    for name,patch in m['endpoint_reconstruction'].items():
        source=OLD/patch['parent_file']
        if sha(source)!=patch['parent_sha256']:raise RuntimeError('Parent changed '+name)
        b=bytearray(source.read_bytes())
        for index,value in patch['byte_replacements']:
            b[index]=value
        (P/'runs'/name).write_bytes(b)
    verify()

def verify():
    m=json.loads((P/'RAW_MANIFEST.json').read_text())
    for name,entry in m['files'].items():
        f=P/'runs'/name
        if f.stat().st_size!=entry['bytes'] or sha(f)!=entry['sha256']:
            raise RuntimeError('Changed raw output '+name)
    print('Verified',len(m['files']),'exact new raw files')

if __name__=='__main__':
    {'pack':pack,'unpack':unpack,'verify':verify}[sys.argv[1]]()

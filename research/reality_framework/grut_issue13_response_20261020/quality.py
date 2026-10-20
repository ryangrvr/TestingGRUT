#!/usr/bin/env python3
"""Response-packet integrity and local-link check; not a physics test."""
from pathlib import Path
import hashlib
import json
import re
import sys

P=Path(__file__).resolve().parent
MAN=P/'PACKET_MANIFEST.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify_links():
    count=0
    for p in P.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in link or link.startswith('#'):continue
            target=(p.parent/link.split('#')[0]).resolve()
            if not target.exists():raise RuntimeError(f'Broken link: {p.name}: {link}')
            count+=1
    return count

if len(sys.argv)>1 and sys.argv[1]=='seal':
    files={f.name:{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(P.iterdir()) if f.is_file() and f!=MAN}
    m={'base_commit':'2ee7f82a62287ae67d92b3b481646f8efed20fc0',
       'diagnostic_freeze':'02712b18bed7519ee4d4f614ac0fd86771a3eddf',
       'review_commit':'55fd11b056385126e876a2b388327eee3028c9d9',
       'files':files,'note':'Manifest excludes itself, extracted runs and compiler products. The raw manifest authenticates all 36 new outputs; original raw hashes are checked by analyze.py. Integrity is not independent scientific replication.'}
    MAN.write_text(json.dumps(m,indent=2)+'\n')
m=json.loads(MAN.read_text())
actual={f.name for f in P.iterdir() if f.is_file() and f!=MAN}
if actual!=set(m['files']):raise RuntimeError('Response file set changed')
for name,item in m['files'].items():
    f=P/name
    if sha(f)!=item['sha256'] or f.stat().st_size!=item['bytes']:raise RuntimeError('Changed file '+name)
print('Integrity:',len(actual),'files;',verify_links(),'local links resolve. Not a scientific replication.')

"""Explicit packaging operation. The read-only verifier never calls this."""
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parents[1]
d={'version':'inflaton-channel-v1','base_commit':'d601fd14009bb7ff99533f575cd047968a65a707',
   'scope':'One-author conditional action/perturbative quantum-wave/averaged receiver calculation; no empirical validation or full stochastic transport',
   'files_sha256':{str(x.relative_to(P)):hashlib.sha256(x.read_bytes()).hexdigest()
       for x in sorted(P.rglob('*')) if x.is_file() and x.name!='MANIFEST.json' and '__pycache__' not in x.parts}}
(P/'MANIFEST.json').write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')

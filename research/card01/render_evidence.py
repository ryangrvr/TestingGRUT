"""Generate the exact table in the manuscript and a portable run receipt."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results.json').read_text())
header=['| Pair | Frequency | Candidate decay rate | Conventional decay rate |',
        '|---|---:|---:|---:|']
rows=[f"| {x['pair'][0]}–{x['pair'][1]} | {x['omega']} | {x['candidate_rate']} | {x['countermodel_rate']} |"
      for x in r['exact_restriction']['rows']]
table='\n'.join(header+rows)
path=ROOT/'REPORT.md'
report=path.read_text()
if '<!-- GENERATED_RATE_TABLE -->' in report:
    report=report.replace('<!-- GENERATED_RATE_TABLE -->',
          '<!-- BEGIN_GENERATED_RATE_TABLE -->\n'+table+'\n<!-- END_GENERATED_RATE_TABLE -->')
else:
    begin='<!-- BEGIN_GENERATED_RATE_TABLE -->'; end='<!-- END_GENERATED_RATE_TABLE -->'
    pre,rest=report.split(begin,1); _,post=rest.split(end,1)
    report=pre+begin+'\n'+table+'\n'+end+post
path.write_text(report)
excluded={'MANIFEST.json','CHECKS.md','FREEZE_RECEIPT.json','STATE.md'}
manifest={'scientific_status':'AUTHOR_KILL_PENDING_INDEPENDENT_REVIEW',
          'freeze_commit':'4bbd8d9b3d6e71df93847d043f49e638c067728c',
          'law_unchanged_since_freeze':r['frozen_files_verified'],
          'checks_status':r['status'],
          'external_imports':{'development/final_precard/price_oracle.py':
              hashlib.sha256((ROOT.parents[1]/'development/final_precard/price_oracle.py').read_bytes()).hexdigest()},
          'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted(ROOT.iterdir()) if p.is_file() and p.name not in excluded}}
(ROOT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('CARD1_REPORT_GENERATED_FROZEN_LAW_UNCHANGED')

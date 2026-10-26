"""Pin inspected inputs and check the additive packet before publication."""
from pathlib import Path
import ast
import hashlib
import json
import re
import subprocess
from audit import ROOT,save

REPO=ROOT.parents[2]
PREFIX='research/reality_framework/'
OLD='225028b56fbc6f11174b8213394874ba6fadad4e'
REVIEW='88fe3b08d890baa0acebcb4f3fb92f487fbe3b10'
F7='f7d641a2570f148b67e4ad05bd8acbea5f97c960'


def main():
    sources=[]
    specs=[]
    for f in ('BACKGROUND_RESULTS.json','BUDGET_RESULTS.json','INPUTS.json','REPORT.md','CLAUDE_HANDOFF.md'):
        specs.append((OLD,PREFIX+'grut_heavy_vector_budget_20261025/'+f,'Read as supplied computed inputs or scope; not rerun.'))
    for p,scope in (
        ('REVIEW_HEAVYWINDOW.md','Full consolidated review read; independent results quoted, not rerun.'),
        ('tracks/hw_critic/CONSOLIDATED.md','Sections 1–8 and verdict tables inspected; code not executed.'),
        ('tracks/hw_gates/code/spotmodel.py','Potential and preparation selectors, approximately lines 1–190; not imported or run.'),
        ('tracks/hw_gates/code/gate1_vector.py','Approximately lines 1–110: spectrum and source definitions, not run. Stale header e-scaling is corrected by the review.'),
    ):
        specs.append((REVIEW,PREFIX+'grut_formation_review_20261025/'+p,scope))
    for p,scope in (
        ('tracks/ls_threshold/code/lsrelic.py','Lines 1–280: state functions, rates, onset and solver structure; not imported or run in this packet.'),
        ('tracks/ls_threshold/code/lsthermo.py','Bose-series thermodynamics read, not imported or run in this packet.'),
    ):
        specs.append((F7,PREFIX+'grut_formation_review_20261025/'+p,scope))
    for commit,path,scope in specs:
        data=subprocess.check_output(['git','show',commit+':'+path],cwd=REPO)
        sources.append(dict(kind='git_file',commit=commit,path=path,sha256=hashlib.sha256(data).hexdigest(),
                            bytes=len(data),read_scope=scope,executed_here=False))
    context=json.loads((ROOT/'ISSUE_CONTEXT.json').read_text())
    for key,comment in context.items():
        sources.append(dict(kind='issue_comment',id=comment['id'],url=comment['url'],
                            raw_body_sha256=hashlib.sha256(comment['body'].encode()).hexdigest(),
                            read_scope='Exact fetched body; context saved in ISSUE_CONTEXT.json.',label=key))
    sources += [
        dict(kind='primary_paper',url='https://arxiv.org/html/2407.12104v2',
             title='Freezing-in Cannibal Dark Sectors',read_scope='Sections 2 and 5.2; moment closure and coupled-species reactions. No model results imported.',executed_here=False),
        dict(kind='primary_paper',url='https://arxiv.org/html/2309.16615v2',
             title='Anatomy of real intermediate state-subtraction scheme',read_scope='Introduction and minimal real-state double-counting discussion. No subtraction scheme implemented here.',executed_here=False)
    ]
    save('SOURCES.json',dict(records=sources,papers_are_context_not_model_validation=True))
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or p.name=='MANIFEST.json' or '__pycache__' in p.parts:continue
        data=p.read_bytes()
        if p.suffix=='.py':ast.parse(data,filename=str(p))
        if p.suffix=='.json':
            json.loads(data,parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
        if p.suffix=='.md':
            for dest in re.findall(r'\[[^\]]*\]\(([^)]+)\)',data.decode()):
                if dest.startswith(('https://','http://','#')):continue
                assert (p.parent/dest.split('#')[0]).exists(),(p,dest)
        files.append(dict(path=str(p.relative_to(ROOT)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    save('MANIFEST.json',dict(files=files,excludes=['MANIFEST.json','__pycache__'],
                              inspected_source_records=len(sources),new_physical_solver=False))
    print('PACKET_OK',len(files),'files',len(sources),'source records')


if __name__=='__main__':main()

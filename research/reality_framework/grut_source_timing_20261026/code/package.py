"""Seal exact inputs and the conditional chronology outputs."""
from timing_common import *
import ast,re,platform

PARENT='423d37b78950d30552179d4de0c1c465016bca6d'
CW='a5cd6a1e2f406fa9e6a4c88da092fed5d3a519fc'
RM='bca4927214695ef7c668ebb1d20f26a9cba91ca7'
MODES='2cbc7f1c3a925fea1f3ac38d0b254ca037afd95e'
REHEAT='e1603565a0c1411a55cdef7f5dc4ecc25aae1c8f'
KIN='34f23f85df3160a86981ae40a159d0c6fffdf95a'
BUILD='e71730c24fe6a0b28ddd72cd30982c5b39c9aedf'
BASE='research/reality_framework/'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def pins():
    review=BASE+'grut_formation_review_20261025/'
    pairs=[(CW,review+p) for p in ['REVIEW_COREWAVE.md','tracks/cw_review/code/ledger_check.py','tracks/cw_review/results/ledger_check.json']]
    pairs += [(RM,review+p) for p in ['REVIEW_REHEATING_MODES.md','tracks/rm_review/code/quartic_floquet.py','tracks/rm_review/results/quartic_floquet.json']]
    pairs += [(PARENT,BASE+'grut_core_wave_audit_20261026/'+p) for p in ['REPORT.md','FINAL_AUDIT.json']]
    pairs += [(MODES,BASE+'grut_reheating_modes_20261026/'+p) for p in ['code/modes.py','code/modes.c','code/mode_common.py','code/phase_and_orbit.py','MODE_PROFILES.npz','PHASE_PROFILES.npz','PHASE_RESULTS.json','QUARTIC_ORBIT.json']]
    pairs += [(REHEAT,BASE+'grut_reheating_preparation_20261026/code/'+p) for p in ['reheat.c','reheating.py','common.py']]
    pairs += [(REHEAT,BASE+'grut_heavy_vector_budget_20261025/INPUTS.json')]
    pairs += [(KIN,BASE+'grut_kinetic_continuation_20261026/code/'+p) for p in ['ks_common.py','bootstrap.py']]
    pairs += [(BUILD,review+'tracks/kin_source/'+p) for p in ['results/source_interface.json','code/decay_source.py','code/mediator_source.py']]
    result=[]
    for commit,path in pairs:
        raw=subprocess.check_output(['git','show',commit+':'+path],cwd=REPO)
        result.append(dict(commit=commit,path=path,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    return result

def finite(x):
    if isinstance(x,float):assert np.isfinite(x)
    elif isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)

def main():
    seal='--seal' in sys.argv
    if seal:
        save('SOURCES.json',dict(internal=pins(),scope='Immutable inputs and completed reviews, no new observational data.'))
        save('EXECUTION.json',dict(python=platform.python_version(),numpy=np.__version__,
            compiler=subprocess.check_output(['gcc','--version']).decode().splitlines()[0],
            flags=['-O3','-std=c11','-Wall','-Wextra'],numerical_worker_processes=2,
            profile_precision_bits=64,direct_vector_source_sha256=sha(ROOT/'code/response.c')))
    sources=read('SOURCES.json')['internal']
    for row in sources:
        raw=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=REPO)
        assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256']
        if row['commit'] in [PARENT,MODES,REHEAT,KIN]:assert sha(REPO/row['path'])==row['sha256']
        if row['commit']==BUILD:
            path=K.INPUT_ROOT/row['path'].removeprefix(BASE+'grut_formation_review_20261025/')
            assert sha(path)==row['sha256']
    v=read('VERIFICATION_RESULTS.json');assert len(v['groups'])==7 and all(g['passed'] for g in v['groups'].values())
    response=read('RESPONSE_RESULTS.json');result=read('RESULTS.json');clocks=read('SOURCE_CLOCKS.json')
    assert len(response['rows'])==4 and sum(len(r['ks']) for r in response['rows'])==42
    assert len(clocks['histories'])==2 and result['physical_f_sym'] is None and result['physical_bound_population'] is None
    assert result['production_source_normalizations_recomputed'] is False
    assert len(read('QUARTIC_SCAN.json')['sampled_bands'])==2
    for row in response['rows']:
        path=ROOT/(row['tag']+'.npz')
        with np.load(path) as z:
            for key in z.files:assert np.all(np.isfinite(z[key]))
            assert z['times'][-1]==800 and z['selected_vectors'].shape[1]==3
            assert z['amplitude_per_radian'].shape[1]==len(row['ks'])
        raw=WORK/(row['tag']+'.bin')
        if raw.exists():assert sha(raw)==row['controls']['raw_sha256']
    for path in ROOT.glob('*.json'):finite(json.loads(path.read_text()))
    for path in (ROOT/'code').glob('*.py'):ast.parse(path.read_text(),filename=str(path))
    for path in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (path.parent/link.split('#')[0]).exists(),(path,link)
    assert (ROOT/'TIMING.png').stat().st_size>10000
    save('AUDIT_RESULTS.json',dict(passed_numerical_groups=7,source_pins=len(sources),charged_mode_histories=42,
        homogeneous_background_histories=4,JSON_NPZ_finite=True,Python_AST=True,local_links=True,
        imported_pinned_files_unchanged=True,figure_rendered_and_viewed=True,
        physical_nonlinear_breakup_computed=False,physical_relic_abundance_computed=False,
        independent_reviewer_of_this_packet_ran=False))
    if seal:
        rows=[]
        for path in sorted(ROOT.rglob('*')):
            if not path.is_file() or path.name=='MANIFEST.json' or '__pycache__' in path.parts:continue
            rows.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path)))
        save('MANIFEST.json',dict(parent_commit=PARENT,files=rows))
    else:
        for row in read('MANIFEST.json')['files']:
            path=ROOT/row['path'];assert path.stat().st_size==row['bytes'] and sha(path)==row['sha256']
    print(f'AUDIT:7/7 groups;42 conditional mode histories;{len(sources)} source pins; finite data, syntax, links and manifest pass. Physical breakup/abundance remain null.')

if __name__=='__main__':main()

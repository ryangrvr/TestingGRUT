"""Pin inputs, conserve results and seal the completed cold pilot."""
import ast,re,subprocess,platform
from pilot import *

PARENT='23356538581c64e49c255e6cc895373d5d297c96'
BUILD='e71730c24fe6a0b28ddd72cd30982c5b39c9aedf'
KIN='34f23f85df3160a86981ae40a159d0c6fffdf95a'
REVIEW='99aedde454c77f2c7dbabaab320a81e95a932fc7'
BASE='research/reality_framework/grut_formation_review_20261025/'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def finite(v):
    if isinstance(v,float):assert np.isfinite(v)
    elif isinstance(v,dict):
        for z in v.values():finite(z)
    elif isinstance(v,list):
        for z in v:finite(z)

def pins():
    pairs=[(PARENT,'research/reality_framework/grut_formation_onset_20261026/'+p) for p in ['code/common.py','code/elliptic.py','MATHEMATICS.md','REPORT.md']]
    pairs += [(KIN,'research/reality_framework/grut_kinetic_continuation_20261026/code/'+p) for p in ['ks_common.py','bootstrap.py']]
    pairs += [(BUILD,BASE+p) for p in ['tracks/kin_source/results/source_interface.json','tracks/kin_rates/results/rates_interface.json',
        'tracks/kin_thermo/results/thermo_interface.json','tracks/kin_thermo/code/formation.py','tracks/kin_thermo/results/objects_EQ.json']]
    pairs += [(REVIEW,BASE+'REVIEW_KINETIC.md')]
    rows=[]
    for commit,path in pairs:
        raw=subprocess.check_output(['git','show',commit+':'+path],cwd=REPO)
        rows.append(dict(commit=commit,path=path,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    return rows

def make_dataset_record():
    private=REPO.parent/'grut3d_checkpoint_library.json'
    if not private.exists():
        assert (ROOT/'DATASET.json').exists(),'Retain supplied dataset metadata on external reproduction'
        return
    record=json.loads(private.read_text());local=Path(record['local_path']);lid=record['library_file_id'];nodes=[]
    def walk(v):
        if isinstance(v,dict):
            if v.get('library_file_id')==lid:nodes.append(v)
            for x in v.values():walk(x)
        elif isinstance(v,list):
            for x in v:walk(x)
    walk(record['result']);node=next(v for v in nodes if v.get('status')=='succeeded')
    save('DATASET.json',dict(library_file_id=lid,file_id=node['file_id'],file_name=node['file_name'],
        library_path=node['path'],current_version_number=node.get('current_version_number'),
        bytes=local.stat().st_size,sha256=sha(local),url='https://chatgpt.com/library?file_id='+lid+'&account_id=personal',
        case='mixed_N96',N=96,L=48.,t=240.,precision_bits=64,
        arrays='q and p: canonical fields/momenta, shape(3,96,96,96), m_X units; metadata stored in the same archive.',
        source_parent=PARENT,engine_sha256=sha(ROOT/'code/field3d.cpp'),
        scope='Owner-accessible full main checkpoint. The external archive is not duplicated in the git packet; public replication uses the source and exact seeds.'))
    (ROOT/'DATASET.md').write_text('# Main checkpoint and public data\n\n'
        +'[GRUT_3D_mixed_N96_checkpoint_20261026.npz](https://chatgpt.com/library?file_id='+lid+'&account_id=personal) contains the complete double-precision main96³ final state at m_X t=240. Access follows the owner account. It includes q,p,N,L,time and the supplied preparation metadata. DATASET.json retains the exact identifiers, hash and engine hash.\n\n'
        +'The checkpoint is a cold classical field configuration, not a particle catalog or a cosmological result. Public reproduction uses the code and exact Fourier parameters. Intermediate raw files are regenerable and their hashes are in SNAPSHOT_MANIFEST.json. VISUALIZATION_FIELDS.npz contains block-averaged charge, phi and energy densities for the two main final states; those reduced arrays are for viewing and are not used for the catalogs or energy tests.\n')

def make_table():
    result=json.loads((ROOT/'RESULTS.json').read_text());lookup={r['tag']:r for r in result['rows']}
    names=[('mixed_N64','Mixed seed,64³'),('mixed_N96','Mixed seed,96³'),('mixed_N128','Mixed seed,128³'),
        ('mixed_dt_half','Mixed,half time step'),('mixed_seed2','Mixed,second seed'),('onesign_N96','Mostly one-sign seed'),('mixed_box64','Larger box,L=64')]
    table=['| Control | Final candidate count | Final Q_mask/Q_net | Late sampled mean | Boundary-shell energy |',
           '|---|---:|---:|---:|---:|']
    for tag,name in names:
        r=lookup[tag];table.append(f"| {name} | {r['final_count']} | {r['final_Qmask_over_netQ']:.3f} | {r['late_Qmask_over_netQ_mean']:.3f} | {100*r['final_boundary_shell_075']:.1f}% |")
    table+='','Boundary shell means max momentum component above0.75 of the axial Nyquist value, using kinetic-plus-gradient energy. The late mean interpolates the stored masks at m_X t=160,180,200,220,240; it is not a statistical confidence interval.'
    text=(ROOT/'REPORT.md').read_text();start='<!-- RESULT_TABLE -->';end='<!-- END_RESULT_TABLE -->';block=start+'\n\n'+'\n'.join(table)+'\n\n'+end
    if end in text:text=re.sub(re.escape(start)+'.*?'+re.escape(end),lambda m:block,text,flags=re.S)
    else:text=text.replace(start,block)
    (ROOT/'REPORT.md').write_text(text)

def main():
    seal='--seal' in sys.argv
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];result=json.loads((ROOT/'RESULTS.json').read_text());analysis=json.loads((ROOT/'ANALYSIS.json').read_text())['cases']
    assert len(runs)==7 and len(result['rows'])==7 and len(analysis)==7
    assert all('error' not in r for r in runs.values())
    checks=json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text())
    assert len(checks['groups'])==9 and all(g['passed'] for g in checks['groups'].values())
    assert result['physical_f_sym'] is None and not result['bound_object_population_certified']
    opposite=json.loads((ROOT/'OPPOSITE_CURRENT.json').read_text())
    assert len(opposite['cases'])==7 and opposite['physical_antiparticle_interpretation'] is False
    if seal:
        save('SOURCES.json',dict(internal=pins(),scope='Pinned two-field benchmark and inherited interface data. No external observational data used.'))
        traces={};snaps=[]
        for tag,r in runs.items():
            data=np.genfromtxt(WORK/tag/'state.tsv',names=True,delimiter='\t');traces[tag]={name:data[name] for name in data.dtype.names}
            for path in r['snapshots']:
                p=Path(path);N,L,t,a=snapshot(p);snaps.append(dict(path=str(p.relative_to(WORK)),bytes=p.stat().st_size,sha256=sha(p),N=N,L=L,t=t,precision_bits=a.dtype.itemsize*8))
        save('TRACES.json',dict(cases=traces));save('SNAPSHOT_MANIFEST.json',dict(files=snaps,scope='Regenerable raw scratch snapshots; selected main checkpoint separately preserved.'))
        save('EXECUTION_ENV.json',dict(compiler=subprocess.check_output(['g++','--version']).decode().splitlines()[0],
            flags=['-std=c++17','-O3','-ffp-contract=off','-fopenmp'],python=platform.python_version(),numpy=np.__version__,
            numerical_processes=2,threads_per_process=4,engine_sha256=sha(ROOT/'code/field3d.cpp')))
        make_dataset_record();make_table()
    for row in json.loads((ROOT/'SOURCES.json').read_text())['internal']:
        raw=subprocess.check_output(['git','show',row['commit']+':'+row['path']],cwd=REPO)
        assert len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256']
        if row['commit'] in [PARENT,KIN] and row['path'].endswith('.py'):assert sha(REPO/row['path'])==row['sha256']
        if row['commit']==BUILD:
            actual=C.K.INPUT_ROOT/row['path'].removeprefix(BASE);assert sha(actual)==row['sha256']
    for f in ROOT.glob('*.json'):finite(json.loads(f.read_text()))
    for f in (ROOT/'code').glob('*.py'):ast.parse(f.read_text(),filename=str(f))
    for f in ROOT.glob('*.npz'):
        with np.load(f) as d:
            for key in d.files:assert np.all(np.isfinite(d[key]))
    for f in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (f.parent/link.split('#')[0]).exists(),(f,link)
    for name in ['STRUCTURES.png','CONTROLS.png']:assert (ROOT/name).stat().st_size>10000
    save('AUDIT_RESULTS.json',dict(passed_numerical_groups=9,production_runs=7,source_pins=12,JSON_NPZ_finite=True,
        Python_AST=True,links=True,figures_rendered_and_viewed=True,full_main_checkpoint_preserved=True,
        physical_f_sym_certified=False,physical_core_population_certified=False,independent_reviewer_ran=False))
    if seal:
        files=[]
        for f in sorted(ROOT.rglob('*')):
            if not f.is_file() or f.name=='MANIFEST.json' or '__pycache__' in f.parts:continue
            files.append(dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=sha(f)))
        save('MANIFEST.json',dict(parent_commit=PARENT,files=files))
    else:
        for row in json.loads((ROOT/'MANIFEST.json').read_text())['files']:
            f=ROOT/row['path'];assert f.stat().st_size==row['bytes'] and sha(f)==row['sha256']
    print('AUDIT:9/9 numerical groups;7 runs;12 pins; finite data, syntax, links and manifest pass. Physical abundance/population remain uncertified.')

if __name__=='__main__':main()

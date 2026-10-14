"""Execute the actual parent and new calculations; audit new outputs.
Does not overwrite recorded outputs except when explicitly invoked with --record.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile,time,platform
import numpy as np
import scipy
P=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    record='--record' in sys.argv
    manifest=json.loads((P/'MANIFEST.json').read_text())
    for name,value in manifest['dependencies'].items():
        assert sha(P/name)==value,('dependency hash',name)
    # Source hashes are immutable; result/verification files are checked separately
    # so an explicit fresh record can replace numerical timing metadata.
    for name,value in manifest['source_sha256'].items():assert sha(P/name)==value,('source hash',name)
    start=time.monotonic();steps=[]
    with tempfile.TemporaryDirectory(prefix='cp12-equal-energy-') as tmp:
        dest=Path(tmp);env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',CP12_OUTPUT_DIR=tmp)
        programs=[('../cp12_preparation_minentropy_audit/check.py','CP12_PREPARATION_MINENTROPY_FORMULA_CHECK_PASS'),
                  ('exact_check.py','CP12_EQUAL_ENERGY_EXACT_UNREGULATED_PASS'),
                  ('preparation.py','CP12_EQUAL_ENERGY_PREPARATION_INITIAL_RESPONSE_PASS'),
                  ('quantum_check.py','CP12_EQUAL_ENERGY_FULL_LAW_PROPAGATION_PASS')]
        for name,marker in programs:
            before=time.monotonic();r=subprocess.run([sys.executable,'-I',str(P/name)],env=env,text=True,capture_output=True)
            assert r.returncode==0,(name,r.stdout[-3000:],r.stderr[-3000:])
            assert marker in r.stdout,(name,'missing pass marker')
            steps.append(dict(script=name,status='PASS',seconds=time.monotonic()-before,stdout_sha256=hashlib.sha256(r.stdout.encode()).hexdigest()))
            print(name,'PASS',flush=True)
        exact=json.loads((dest/'exact_results.json').read_text())
        assert [F(r['full_coefficients'][2]) for r in exact['rows']]==[-F(1,810),F(1,1215)]
        assert {r['mean_energy'] for r in exact['rows']}=={'32/15'}
        assert {r['variance'] for r in exact['rows']}=={'5402/675'}
        big=json.loads((dest/'preparation_results.json').read_text())
        assert len(big['rows'])==18 and all(r['Q']==3*r['N'] for r in big['rows'])
        last={r['pattern']:r for r in big['rows'] if r['N']==256}
        assert last['concentrated']['observables']['squares']['T2_coefficient']>0
        assert last['concentrated']['observables']['wedges']['T2_coefficient']>0
        assert last['spread']['observables']['squares']['T2_coefficient']<0
        quantum=json.loads((dest/'quantum_results.json').read_text())
        last_comp=quantum['cutoff_comparisons'][-1]
        assert last_comp['cutoffs']==[4,5]
        differences=[d for row in last_comp['rows'] for d in row['absolute_differences'].values()]
        assert max(d['squares'] for d in differences)<2e-9
        assert max(d['electric_energy'] for d in differences)<2e-5
        small=json.loads((dest/'quantum_L5.json').read_text())
        assert small['dimension']==10275
        final=small['finite_time'][-1]['preparations']
        assert final['spread']['squares']>.25 and final['concentrated']['squares']<.12
        assert .20<final['dephased']['squares']<.21
        source=small['pure_initial_coefficients']
        for row in exact['rows']:
            c=source[row['pattern']]['observables']['squares']['full_T_coefficients']
            assert max(abs(float(F(x))-c[j]) for j,x in enumerate(row['full_coefficients']))<1e-11
        regenerated={f.name:sha(f) for f in sorted(dest.glob('*.json'))}
        if record:
            for f in dest.glob('*.json'):(P/f.name).write_bytes(f.read_bytes())
        out=dict(status='CP12_ENERGY_MATCHED_PREPARATION_VERIFIED',python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,steps=steps,
                 checks=['Pinned parent checker executed verbatim','Exact unregulated full-H certificates','Independent rational and quotient high-density derivatives','Full coupled pure/mixed evolution at five cutoffs','Observable-specific cutoff convergence','Exact-versus-matrix coefficient agreement'],
                 seconds=time.monotonic()-start,regenerated_sha256=regenerated,
                 source_sha256={name:sha(P/name) for name in manifest['source_sha256']},
                 limits=['No full Q=3N finite-time propagation','No infinite-cutoff finite-time enclosure','No controller dynamics or state selection','Square counts are not disjoint packing or spatial volume growth'])
    if record:(P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()

"""Fresh execution of pinned parent and all four cubic response controls."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,platform,subprocess,sys,tempfile,time
import numpy as np
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    manifest=json.loads((P/'MANIFEST.json').read_text())
    for name,h in manifest['dependencies'].items():assert sha(P/name)==h,('dependency',name)
    for name,h in manifest['source_sha256'].items():assert sha(P/name)==h,('source',name)
    started=time.monotonic();steps=[];record='--record' in sys.argv
    with tempfile.TemporaryDirectory(prefix='cp12-cubic-topology-') as tmp:
        folder=Path(tmp);env=dict(os.environ,CP12_OUTPUT_DIR=tmp,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        programs=[('../cp12_phase_motif_indistinguishability/check.py','CP12_PHASE_MOTIF_TOPOLOGY_BLINDNESS_EXACT_PASS'),
                  ('symbolic_check.py','CP12_TOPOLOGY_CUBIC_SYMBOLIC_PASS'),
                  ('full_moments.py','CP12_TOPOLOGY_FULL_H_UNREGULATED_CUBIC_PASS'),
                  ('topology.py','CP12_DEGREE_MATCHED_TOPOLOGY_DERIVATIVES_PASS'),
                  ('asymptotic_check.py','CP12_TOPOLOGY_ASYMPTOTIC_FORMULA_PASS')]
        for name,marker in programs:
            begin=time.monotonic();r=subprocess.run([sys.executable,'-I',str(P/name)],env=env,text=True,capture_output=True)
            assert r.returncode==0,(name,r.stdout[-2000:],r.stderr[-2000:])
            assert marker in r.stdout,(name,'missing pass marker')
            steps.append(dict(script=name,status='PASS',seconds=time.monotonic()-begin,stdout_sha256=hashlib.sha256(r.stdout.encode()).hexdigest()))
            print(name,'PASS',flush=True)
        small=json.loads((folder/'full_moments_results.json').read_text())
        A,B=small['rows'];assert small['Q']==12
        assert A['mean_energy']==B['mean_energy']=='156/25'
        assert A['variance']==B['variance']=='780063/32500'
        for obs in ['triangles','squares','wedges']:assert A['full_coefficients'][obs][:3]==B['full_coefficients'][obs][:3]
        assert F(A['full_coefficients']['squares'][3])==F(1243,273000000)
        assert F(B['full_coefficients']['squares'][3])==-F(37,9100000)
        assert A['full_coefficients']['wedges']==B['full_coefficients']['wedges']
        assert all(F(r['noninvariant_cubic_control']['full_minus_neutral'])==-F(4,375) for r in [A,B])
        assert all(r['phase_conjugation']['real_H_action_verified'] for r in [A,B])
        big=json.loads((folder/'topology_results.json').read_text())
        assert [r['N'] for r in big['models']]==[24,48,96,192]
        assert all(r['Q']==3*r['N'] and r['observables']['wedges']['B_minus_A_T3_over_sin']=='0' for r in big['models'])
        assert all(r['observables']['triangles']['B_minus_A_T3']>0 and r['observables']['squares']['B_minus_A_T3']<0 for r in big['models'])
        symbolic=json.loads((folder/'symbolic_results.json').read_text());assert symbolic['probability_polynomial_degree']==3
        asymptotic=json.loads((folder/'asymptotic_results.json').read_text())
        for b,a in zip(big['models'],asymptotic['rows']):
            assert b['N']==a['N']
            assert abs(b['observables']['triangles']['B_minus_A_T3']-a['triangle_difference_T3'])<1e-18
            assert abs(b['observables']['squares']['B_minus_A_T3']-a['square_difference_T3'])<1e-18
        outputs={f.name:sha(f) for f in sorted(folder.glob('*.json'))}
        if record:
            for f in folder.glob('*.json'):(P/f.name).write_bytes(f.read_bytes())
    out=dict(status='CP12_CUBIC_TOPOLOGY_VERIFIED',python=platform.python_version(),numpy=np.__version__,steps=steps,
             checks=['Byte-identical parent checker run','Formal cubic polynomial identity','Independent full unregulated Hamiltonian moments',
                     'Non-invariant hopping correction','Matched degrees and energy moments','Independent factorial and overlap formulas',
                     'Same Q=3N sector','Controlled N^-5 and N^-6 derivative asymptotics','Full-H conjugation covariance'],
             seconds=time.monotonic()-started,regenerated_sha256=outputs,source_sha256=manifest['source_sha256'],
             limits=['No high-density finite-time dynamics','No uniform Taylor remainder','No disjoint-loop packing, dimension or propagation test',
                     'No autonomous preparation or state-selection mechanism','Auxiliary Q=2N is not main Q=3N'])
    if record:(P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()

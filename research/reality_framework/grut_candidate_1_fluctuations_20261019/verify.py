#!/usr/bin/env python3
"""Check linear algebra, mode controls, pole calculations and scope records."""
from pathlib import Path
import importlib.util
import hashlib,json,csv
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.sparse import diags

HERE=Path(__file__).resolve().parent
REF=HERE.parent/'grut_candidate_1_20261018/code/results_profile.npz'


def main():
    checks=[]
    def check(name,value,detail=None):checks.append({'name':name,'passed':bool(value),'detail':detail})
    m=json.loads((HERE/'mode_results.json').read_text())
    p=json.loads((HERE/'pole_results.json').read_text())
    n=json.loads((HERE/'neutral_pole_results.json').read_text())
    a=json.loads((HERE/'angular_bound_results.json').read_text())
    check('original profile SHA256 unchanged',hashlib.sha256(REF.read_bytes()).hexdigest()=='0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3')
    check('all40 prescribed spectra completed',m['completed'] and len(m['spectra'])==40 and not m['skipped'])
    check('no linear eigensolver exception or nonconvergence',all(t.get('status','converged')=='converged' for t in m['attempts']))
    check('all discrete stationary solves accepted',all(t.get('success',True) for t in m['attempts']))
    check('vacuum charged and neutral dispersions',max(v['maximum_charged_eta_frequency_difference'] for v in m['vacuum_checks'])<1e-10)
    check('vacuum D dispersions',max(v['maximum_D_frequency_difference'] for v in m['vacuum_checks'])<1e-10)
    check('no resolved nonsymmetry growth at declared accuracy',max(s['max_real_part_excluding_symmetry'] for s in m['spectra'])<1e-8)
    check('returned eigenvector residuals controlled',max(s['max_generator_residual'] for s in m['spectra'])<1e-8)
    check('nonzero returned modes satisfy linear charge',max(s['max_deltaQ_normalized_nonzero'] for s in m['spectra'])<1e-9)
    check('no nonsymmetry internal mode in declared finite search window',all(not s['internal_positive_frequencies'] for s in m['spectra']))
    check('radial fixed-charge amplitude energies positive',min(s['energy_checks']['real_amplitude_lowest_energy_eigenvalues'][0] for s in m['spectra'] if s['ell']==0)>.35)
    check('ell2 amplitude energies positive',min(s['energy_checks']['real_amplitude_lowest_energy_eigenvalues'][0] for s in m['spectra'] if s['ell']==2)>.25)
    check('phase energies nonnegative within zero roundoff',min(s['energy_checks']['phase_lowest_energy_eigenvalues'][0] for s in m['spectra'])>-1e-9)
    check('D energies bounded away from zero',min(s['energy_checks']['D_lowest_frequency'][0] for s in m['spectra'])>=1)
    tight=next(s for s in m['spectra'] if s['method']=='dense_phase_tight')
    for key in ['phase_generator_residual','charge_left_null_residual','phase_Jordan_residual']:
        check('tightened '+key,tight['energy_checks'][key]<1e-10)
    check('phase Jordan partner is charge changing',abs(tight['energy_checks']['Jordan_partner_deltaQbar']+tight['energy_checks']['dQbar_dOmega_discrete'])<1e-8)
    initial=json.loads((HERE/'mode_initial_results.json').read_text())
    check('initial apparent phase growth retained',any(any(root['s_real']>1e-7 and root['symmetry']=='phase' for root in s['positive_real_roots']) for s in initial['spectra']))
    trans=[]
    for N in [120,180,400,800]:
        s=next(x for x in m['spectra'] if x['N']==N and x['L']==30 and x['ell']==1 and x['method'] in ['dense','sparse'])
        roots=[r for r in s['smallest_roots'] if r['symmetry']=='translation' and r['s_imag']>0]
        trans.append(min(r['s_imag'] for r in roots))
    check('translation frequencies approach zero with refinement',all(trans[i]>trans[i+1] for i in range(3)),trans)
    def matched(ell,L):
        return next(s for s in m['spectra'] if s['ell']==ell and s['L']==L and s['method'] in ['matched_box','matched_ell5'])
    for ell in range(6):
        s30=matched(ell,30);s40=matched(ell,40)
        def lowest(s):return min(r['s_imag'] for r in s['smallest_roots'] if r['s_imag']>1e-5 and r['symmetry'] is None)
        check('ell'+str(ell)+' box level approaches charged continuum',.2<lowest(s40)<lowest(s30))
    t30=min(r['s_imag'] for r in matched(1,30)['smallest_roots'] if r['symmetry']=='translation' and r['s_imag']>0)
    t40=min(r['s_imag'] for r in matched(1,40)['smallest_roots'] if r['symmetry']=='translation' and r['s_imag']>0)
    check('translation matched-h volume sensitivity small',abs(t30-t40)<1e-8)
    check('radial grand-canonical negative direction preserved',matched(0,40)['energy_checks']['unconstrained_amplitude_eigenvalues'][0]<0)
    check('all10 X pole checks passed',p['all_checks_pass'] and len(p['checks'])==10)
    check('all6 tensor pole checks passed',n['all_checks_pass'] and len(n['checks'])==6)
    for filename,result in [('pole_protocol.json',p),('neutral_pole_protocol.json',n),('angular_bound_protocol.json',a)]:
        check(filename+' matches executed input',hashlib.sha256((HERE/filename).read_bytes()).hexdigest()==result['protocol_sha256'])
    check('angular comparison explicitly not half-line certified',a['half_line_bound_certified'] is False)
    check('conditional angular cutoff consistent',a['conditional_first_excluded_ell']==6 and 6*7>a['numerical_margin_comparison_value'])
    # Independent canonical/velocity-generator algebra on a small sampled
    # background. No stationary solver or frozen execution output is modified.
    spec=importlib.util.spec_from_file_location('grut_mode_operators',HERE/'modes.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    original=np.load(REF);spline=CubicSpline(original['rho'],original['y'],axis=1);O=float(original['Omega'])
    N=24;L=20.;h=L/(N+1);r=h*np.arange(1,N+1);F,_,S,_=spline(r)
    lap=diags([-np.ones(N-1)/h**2,2*np.ones(N)/h**2,-np.ones(N-1)/h**2],[-1,0,1],format='csr')
    bg={'N':N,'L':L,'h':h,'r':r,'F':F,'S':S,'D':lap}
    A,K,G,Lu,Lw,Le,LD,C=mod.operators(bg,0,O);A=A.toarray();K=K.toarray();G=G.toarray();I=np.eye(3*N);Z=np.zeros_like(I)
    W=np.block([[K,Z],[Z,I]])
    check('gyroscopic G skew symmetric',np.max(abs(G+G.T))<1e-14)
    check('stiffness K symmetric',np.max(abs(K-K.T))<1e-14)
    check('quadratic energy conserved by full generator',np.max(abs(A.T@W+W@A))<1e-12)
    J=G/(2*O);PX=-J@J;HV=K+O*O*PX
    T=np.block([[I,Z],[O*J,I]]);Ti=np.block([[I,Z],[-O*J,I]])
    canonical=np.block([[-O*J,I],[-HV,-O*J]])
    check('canonical and velocity generators similar',np.max(abs(T@A@Ti-canonical))<1e-12)
    pr=r*F;c=np.r_[2*O*pr,np.zeros(3*N),-pr,np.zeros(N)]
    expected=np.r_[np.zeros(N),Lw.toarray()@pr,np.zeros(4*N)]
    check('charge conservation defect equals background phase defect',np.max(abs(c@A-expected))<1e-12)
    K2=mod.operators(bg,2,O)[1].toarray();K8=mod.operators(bg,8,O)[1].toarray()
    check('higher-ell stiffness centrifugal identity',np.max(abs(K8-K2-np.diag(np.tile((8*9-2*3)/r**2,3))))<1e-12)
    with (HERE/'mode_eigenvalues.csv').open() as f:rows=list(csv.DictReader(f))
    check('raw finite-grid eigenvalues preserved',len(rows)==sum(s['eigenvalues_returned'] for s in m['spectra']))
    result={'checks_passed':sum(c['passed'] for c in checks),'total_checks':len(checks),'all_checks_pass':all(c['passed'] for c in checks),'checks':checks,
            'limits':'Internal algebra/controls only; no certified continuum spectrum, renormalized finite-Q quantum energy or lifetime.'}
    (HERE/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['checks_passed','total_checks','all_checks_pass']}))
    for c in checks:
        if not c['passed']:print(json.dumps(c))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()

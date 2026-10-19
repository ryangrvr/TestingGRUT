#!/usr/bin/env python3
"""Recompute algebra/normalization and audit stored numerical evidence."""
import contextlib
import hashlib
import io
import json
import shutil
import tempfile
from pathlib import Path
import numpy as np
from scipy.integrate import quad
import graphs
import spectral as sp
import theory_toy

HERE=Path(__file__).resolve().parent


def load(name):return json.loads((HERE/name).read_text())


def main():
    checks=[]
    def check(name, passed, detail=None):checks.append({'name':name,'passed':bool(passed),'detail':detail})
    g=load('graphs_results.json'); s=load('spectral_results.json'); a=load('audit_results.json'); t=load('tail_refinement_results.json'); b=load('binding_results.json'); toy=load('theory_toy_results.json')
    check('original profile hash',hashlib.sha256(sp.REF.read_bytes()).hexdigest()=='0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3')
    for name,data in [('graphs',g),('audit',a),('tail_refinement',t),('theory',toy)]:
        check(name+' protocol identity',hashlib.sha256((HERE/(name+'_protocol.json')).read_bytes()).hexdigest()==data['protocol_sha256'])
    check('toy source identity',hashlib.sha256((HERE/'theory_toy.py').read_bytes()).hexdigest()==toy['implementation_sha256'])
    for name,data in [('graphs',g),('audit',a),('tail_refinement',t),('binding',b)]:check(name+' controls passed',data['all_checks_pass'])
    check('complete original angular-frequency grids',s['complete_n6'] and s['complete_n8'] and not s['skipped'] and len(s['partial_integrals'])==82)
    check('all original numerical points succeeded',all(x['success'] for x in s['points']))
    check('complete stricter replay',a['completed'] and len(a['replay_points'])==1312 and len(a['replay_partial_integrals'])==41)
    check('all sampled determinant phases small',max(abs(x['logJost_phase_mod_2pi']) for x in s['points']+a['replay_points'])<1e-7)
    check('all sampled Born imaginary parts small',max(abs(x[k]) for x in s['points']+a['replay_points'] for k in ['Born1_imag','Born2_imag'])<1e-10)
    check('quadrature 6/8 comparison',abs(s['energy_n6']-s['energy_n8'])<1e-4)
    check('radial-domain-origin replay within declared allowance',abs(a['replay_difference'])<.02)
    check('replay summation arithmetic',abs(sum(x['energy'] for x in a['replay_partial_integrals'])-a['energy_replay'])<1e-12)
    for ell,zeta in [(0,.3),(2,.7)]:
        rows=sorted([x for x in a['independent_fd'] if x['ell']==ell and x['zeta']==zeta],key=lambda x:x['n'])
        errors=[abs(x['discrepancy']) for x in rows]
        check('independent FD convergence ell'+str(ell),len(rows)==3 and errors[0]>3*errors[1]>9*errors[2] and errors[-1]<.002,errors)
    # A genuinely fresh finite canonical integral/exact-polynomial run, in temp.
    with tempfile.TemporaryDirectory(prefix='grut-energy-toy-') as tmp:
        old=theory_toy.HERE; theory_toy.HERE=Path(tmp);shutil.copy2(HERE/'theory_protocol.json',Path(tmp)/'theory_protocol.json')
        try:
            with contextlib.redirect_stdout(io.StringIO()):theory_toy.main()
            fresh=json.loads((Path(tmp)/'theory_toy_results.json').read_text())
        finally:theory_toy.HERE=old
    check('fresh exact fixed-charge polynomial',fresh['exact_full_minus_zeta2_reduced_residual']=={})
    check('fresh toy full-vacuum logdet/frequency normalization',abs(fresh['integral_frequency_difference'])<2e-9)
    check('fresh toy reproduces stored answer',abs(fresh['logdet_integral_over_2pi']-toy['logdet_integral_over_2pi'])<1e-12)
    check('wrong phase-coordinate deletion is detected',fresh['wrong_reduced_vs_full_vacuum_logdet_UV_power']==-2)
    check('unfixed-charge negative curvature retained',fresh['two_coordinate_toy']['unfixed_charge_curvature_exact'].startswith('-'))
    for ell,zeta in [(0,.3),(2,.7),(40,5.)]:
        v=sp.run_point(ell,zeta,scale=0.,R=25.,r0=.001,tol=1e-10,label='fresh_verifier_vacuum')
        check('fresh continuum vacuum ell'+str(ell),v['success'] and abs(v['total_remainder'])<1e-7,v['total_remainder'])
    # Independent adaptive Feynman-parameter quadrature, not the stored graph grid.
    specs=[(2,2,0),(1,1,0),(1,1,4*sp.OMEGA**2),(2,1,sp.OMEGA**2),(1.01,1.01,0)]
    kk=np.array([0.,.3,1.,3.]); ker,_=graphs.kernels(kk,sp.OMEGA,128)
    error=0.
    for i,(ma,mb,freq) in enumerate(specs):
        for j,k in enumerate(kk):
            exact=quad(lambda x:-np.log(x*ma+(1-x)*mb-x*(1-x)*(freq-k*k)),0,1,epsabs=1e-12)[0]
            error=max(error,abs(exact-ker[j,i]))
    check('independent adaptive bubble-kernel check',error<1e-11,error)
    for R in [20.,25.]:
        restored=next(x for x in a['support_graphs'] if x['R']==R)
        check('graph support agreement R'+str(R),abs(restored['first']+restored['second']-g['selected_finest']['restored_graph_total'])<1e-6)
    budget=b['error_budget_over_gv']
    check('tail and empirical allowances all included',budget['rounded_reported_allowance']>=budget['sum_before_roundup']>budget['frequency_tail_analytical']+budget['angular_tail_analytical'])
    check('integer target not floating charge',b['N_integer']==39944 and abs(b['Q_reference_continuous']-b['N_integer'])>.1)
    check('integer state solved independently',b['integer_check']['integer_BVP_success'] and abs(b['integer_check']['Qbar_at_N']/.01-b['N_integer'])<1e-5)
    check('integer energy linearization checked',abs(b['direct_integer_minus_linear_over_gv'])<1e-5)
    check('binding arithmetic',abs(b['soliton_0plus1_energy_over_gv']-b['free_NX_threshold_over_gv']-b['DeltaE_freeX_over_gv'])<1e-10)
    check('free-carrier range strictly negative',b['DeltaE_freeX_numerical_range'][1]<0)
    check('reported range explicitly qualified','Not a certified exact-profile' in budget['qualification'])
    result={'checks_passed':sum(x['passed'] for x in checks),'total_checks':len(checks),'all_checks_pass':all(x['passed'] for x in checks),'checks':checks,'scope':'Internal checks and fresh toy/vacuum/kernel calculations. Full stricter continuum replay is in audit_results.json. Not independent external validation, exact-profile interval certification or a lifetime calculation.'}
    (HERE/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()

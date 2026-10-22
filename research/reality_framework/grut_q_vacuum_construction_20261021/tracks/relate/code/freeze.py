"""Freeze a scoped, versioned comparison; never marks exposed data as held out.

Offline: reads recorded SMDR results and computes covariance-aware replay.
Optional --smdr-home: uses freshly built SMDR through the pinned C driver;
no imports from the previous Python matching or RG implementations. Fits a
matching interpolant to recorded couplings, then roots beta_lambda=0 along
each flow (rather than a parabolic minimum on the original 241-point grid).
Also checks one direct full SMDR fit at the resulting critical point.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import subprocess

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parents[1]
ROOT=HERE.parents[1]
I=json.loads((HERE/'INPUTS.json').read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def smdr_check(home):
    home=Path(home).resolve()
    grid=json.loads((ROOT/'tracks/v_vacua/results/v3a_smdr_fit_grid.json').read_text())['fits']
    names=['g3','g','gp','yt','yb','ytau','k']
    b=[x for x in grid if abs(x['inputs'][1]-125.2)<1e-10 and abs(x['inputs'][2]-.118)<1e-10]
    co={k:np.polyfit([x['inputs'][0] for x in b],[x[k] for x in b],2) for k in names}
    def get(mh,als):
        return next(x for x in grid if np.allclose(x['inputs'],[171.2,mh,als],rtol=0,atol=1e-10))
    dm={k:get(125.7,.118)[k]-get(124.7,.118)[k] for k in names}
    da={k:(get(125.2,.1189)[k]-get(125.2,.1171)[k])/.0018 for k in names}

    def matching(mt,mh,als):
        return {k:float(np.polyval(co[k],mt)+dm[k]*(mh-125.2)+da[k]*(als-.118)) for k in names}

    def call(args):
        p=subprocess.run([str(home/'smdr_driver')]+[str(x) for x in args],cwd=home,
                         check=True,capture_output=True,text=True)
        return p.stdout

    def run(mt,c,order=3,n=801):
        txt=call(['run',order,format(mt,'.17g')]+[format(c[k],'.17g') for k in names]+[n,1e8,1.2209e19])
        rows=[]
        for line in txt.splitlines():
            if line.startswith('pt '):
                w=line.split();rows.append({w[j]:float(w[j+1]) for j in range(1,len(w),2)})
        x=np.log([r['Q'] for r in rows]);kl=CubicSpline(x,[r['k'] for r in rows])
        bl=CubicSpline(x,[r['betak'] for r in rows])
        opts=[float(x[0]),float(x[-1])]
        for j in range(len(rows)-1):
            if rows[j]['betak']<0<rows[j+1]['betak']:
                opts.append(brentq(bl,float(x[j]),float(x[j+1]),xtol=1e-12))
        xx=min(opts,key=lambda a:float(kl(a)))
        return float(kl(xx)),math.exp(xx)

    def critical(mh=125.2,als=.118,order=3):
        mt=brentq(lambda m:run(m,matching(m,mh,als),order)[0],170,172.5,xtol=2e-8)
        return [mt,run(mt,matching(mt,mh,als),order)[1]]
    vals={str(k):critical(order=k) for k in [3,4,5]}
    slopes={'Mh':critical(mh=125.7)[0]-critical(mh=124.7)[0],
            'alpha_s_per_0.0009':(critical(als=.1189)[0]-critical(als=.1171)[0])/2}
    mt=vals['3'][0]
    text=call(['fit',format(mt,'.17g'),125.2,.118]);cdir={}
    for line in text.splitlines():
        if line.startswith('fit '):
            w=line.split();cdir={w[j]:float(w[j+1]) for j in range(1,len(w),2)}
    if not cdir: raise RuntimeError('No full SMDR fit output')
    kr,qr=run(mt,cdir)
    lo,hi=mt-.01,mt+.01
    slope=(run(hi,matching(hi,125.2,.118))[0]-run(lo,matching(lo,125.2,.118))[0])/.02
    return {'critical_by_loop_order':vals,'slopes':slopes,
            'minimum_from_direct_matching_at_critical_mass':kr,'direct_matching_minimum_scale_GeV':qr,
            'direct_matching_mass_correction_estimate_GeV':-kr/slope,
            'driver_sha256':sha(ROOT/'tracks/v_vacua/code/smdr_driver.c'),
            'reference_model_sha256':sha(home/'ReferenceModel.dat'),
            'library_sha256':sha(home/'libsmdr.a'),
            'scope':'Reuses published SMDR library and archived fit couplings; new minimization and direct matching check, not an independent beta-function calculation.'}


def evaluate(mt,mh,als,C,coeff,extra_sigma=0):
    ref=I['fit_reference']
    crit=coeff['intercept_GeV']+coeff['dMt_dMh']*(mh-ref['Mh_GeV'])+coeff['dMt_per_alpha_s_unit']*(als-ref['alpha_s_5_MZ'])/ref['alpha_s_unit']
    # Full covariance order is (Mt, Mh, alpha_s). Cross terms remain available.
    grad=np.array([1,-coeff['dMt_dMh'],-coeff['dMt_per_alpha_s_unit']/ref['alpha_s_unit']])
    em=I['illustrative_error_model']
    variance=float(grad@np.asarray(C)@grad)+em['matching_sigma_GeV']**2+em['criterion_sigma_GeV']**2+extra_sigma**2
    z=(mt-crit)/math.sqrt(variance)
    return {'critical_top_mass_GeV':crit,'difference_GeV':mt-crit,'sigma_GeV':math.sqrt(variance),
            'indicative_z':z,'Gaussian_one_sided_tail_if_error_model_adopted':math.erfc(z/math.sqrt(2))/2,
            'formal_significance':False}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--smdr-home');args=ap.parse_args()
    p=ROOT/'tracks/v_vacua/results/v3b_sm_check.json'
    raw=json.loads(p.read_text())['C_smdr_matching_smdr_running']
    co={'intercept_GeV':raw['Mtcrit_PDG2025']['loopOrder3'][0],
        'dMt_dMh':raw['dMtcrit/dMh'],'dMt_per_alpha_s_unit':raw['dMtcrit per 0.0009 as']}
    obs=I['reference_input_snapshot'];ref=I['fit_reference']
    C=np.diag([obs['direct_mass_sigma_GeV']**2,obs['Mh_sigma_GeV']**2,obs['alpha_s_sigma']**2])
    direct=evaluate(obs['direct_mass_GeV'],ref['Mh_GeV'],ref['alpha_s_5_MZ'],C,co)
    mc=evaluate(obs['direct_mass_GeV'],ref['Mh_GeV'],ref['alpha_s_5_MZ'],C,co,.3)
    Cx=C.copy();Cx[0,0]=obs['xsec_pole_mass_sigma_GeV']**2
    xsec=evaluate(obs['xsec_pole_mass_GeV'],ref['Mh_GeV'],ref['alpha_s_5_MZ'],Cx,co)
    checks={'replay_direct_z':abs(direct['indicative_z']-3.139089421156063)<1e-10,
            'replay_MC_sensitivity_z':abs(mc['indicative_z']-2.675002597760182)<1e-10,
            'replay_xsec_z':abs(xsec['indicative_z']-1.7270296842279105)<1e-10}
    result={'version':I['version'],'recorded_utc_date':I['recorded_utc_date'],
        'status':'Conditional relation recorded after exposure to all checkpoint data; not a new or GRUT-specific prediction',
        'conditions':['One shared mu; minimally coupled, q-independent matter; chosen stable q branch',
            'P_q matched to a specified true final vacuum, rather than independently to each phase',
            'SM Higgs potential describes all relevant competing vacua up to the stated Planck cutoff',
            'No relevant new thresholds or unsupplied higher-dimensional operators; conventional flat-space SM stability proxy',
            'Small present vacuum gap is separately demanded; the action does not select that gap'],
        'branches':{
            'EW_is_true_final_vacuum':{'relation':'Mt <= Mt_critical, within the declared theoretical approximation',
                 'vacuum_density':'zero plus transient corrections; retains the failure to produce observed positive dark energy'},
            'EW_is_false_near_degenerate_vacuum':{'relation':'V_EW - V_true = rho_Lambda,obs; near-critical equality rather than absolute stability',
                 'gap':'supplied condition; no mechanism derives its magnitude',
                 'precision':'Finite-order SM running cannot resolve a ~1e-120 quartic splitting'},
            'MPP_coexistence':{'relation':'lambda=beta_lambda=0 at a high scale',
                 'status':'Additional Froggatt-Nielsen assumption; not implied by q dynamics'}},
        'operational_RG_proxy':I['nominal_criterion'],'local_linearization':co,
        'local_use_box':{'Mh_GeV':[124.7,125.7],'alpha_s_5_MZ':[.1171,.1189],
                         'outside_box':'Run the recorded RG criterion, not an extrapolation of this line.'},
        'higher_order_intercepts_GeV':{k:raw['Mtcrit_PDG2025']['loopOrder'+k][0] for k in ['3','4','5']},
        'retrospective_snapshot_arithmetic':{'direct_mass_assumed_pole':direct,'direct_with_MC_to_pole_sensitivity':mc,'cross_section_pole':xsec},
        'future_use':{'measurement':'First subsequently published qualifying e+e- top-threshold joint mass and alpha_s likelihood, with Higgs input and theory/mass-scheme treatment stated',
            'eligible_data_rule':'Only data unobserved by this analysis when this version was published; no retrospective reclassification',
            'primary_mass':'Short-distance mass consistently matched to MSbar couplings; never silently identify a Monte Carlo mass with a pole mass',
            'covariance':'Use the full joint covariance or likelihood; gradient is (1,-dMt_dMh,-dMt_dalpha_s)',
            'decision':'A discrepancy concerns the named conditional SM/final-vacuum hypothesis. A 5-sigma flag under the illustrative Gaussian model requires a separately justified likelihood; no fixed future significance is forecast.',
            'support':'Agreement also supports ordinary SM stability and cannot uniquely support GRUT',
            'revision_rule':'New physics, matching or criterion changes require a new version; retain this version and report both.'},
        'provenance':{'source_commit':I['source_commit'],'source_results_sha256':sha(p),'inputs_sha256':sha(HERE/'INPUTS.json'),'freeze_code_sha256':sha(Path(__file__))}}
    if args.smdr_home:
        result['SMDR_reproduction']=smdr_check(args.smdr_home)
        rep=result['SMDR_reproduction']
        checks['new_minimizer_reproduces_critical_mass']=abs(rep['critical_by_loop_order']['3'][0]-co['intercept_GeV'])<.001
        checks['direct_full_matching_agrees_near_boundary']=abs(rep['direct_matching_mass_correction_estimate_GeV'])<.001
    else:
        old=HERE/'results/frozen_predictions.json'
        if old.exists():
            prev=json.loads(old.read_text())
            if prev.get('provenance')==result['provenance'] and 'SMDR_reproduction' in prev:
                result['SMDR_reproduction']=prev['SMDR_reproduction']
                result['SMDR_reproduction_note']='Retained prior computed record; library was not re-executed by this offline replay.'
    result['checks']=checks
    result['checks_passed']=all(checks.values())
    if not result['checks_passed']:raise AssertionError(checks)
    (HERE/'results/frozen_predictions.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'checks':checks,'conditional_critical_mass_GeV':co['intercept_GeV'],
                      'indicative_z':[direct['indicative_z'],mc['indicative_z'],xsec['indicative_z']]},indent=2))


if __name__=='__main__':main()

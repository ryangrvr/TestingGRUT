#!/usr/bin/env python3
"""Deciding analytic, profile, accounting and detector numerical controls.

One author. Reimplementation of formulas is not independent scientific review.
The optional public-toy checks rerun the likelihood at a refined strength grid.
"""
import argparse
import json

import numpy as np
from scipy.integrate import simpson

from widths import ROOT, Widths, HBAR, hidden, visible_coefficient, population
from scan import Object, detector, row_for, PUB, base
from fixed_yield import load_toys, original, YQ_OLD, BDM


def run(cache=None):
    checks = []
    diagnostics = {}
    def check(name, condition, evidence=None):
        checks.append(dict(name=name, passed=bool(condition), evidence=evidence))

    # Constant-amplitude invariant two-body phase-space, distinct/identical.
    differences = []
    for v in (.5, 1., 7.):
        for lam in (.005, .02, .0202, .021, .03, .04):
            ms = np.sqrt(2*lam)*v
            expected = []
            for mass, vertex, symmetry in ((.1*v, .02*v, 1),
                                           (np.sqrt(.0101)*v, .0002*v, 2)):
                beta = np.sqrt(max((ms-2*mass)*(ms+2*mass), 0.))/ms
                if abs(ms-2*mass) < 1e-14*ms:
                    beta = 0.
                expected.append(vertex**2/(2*ms)*beta/(8*np.pi)/symmetry)
            calc = hidden(v, lam)
            differences.extend(abs(c-e)/max(e, 1e-100) for c,e in zip(calc,expected))
    check('Hidden widths equal invariant two-body phase-space including identical DD',
          max(differences) < 2e-12, max(differences))
    check('Vacuum hidden-decay and homogeneous-binding intervals overlap',
          2*.1**2 < .1**4/.001 and hidden(1,.02)[0] == 0 and hidden(1,.0202)[1] == 0)

    widths = Widths()
    me, vh, ms = .00051099895, 246., np.sqrt(.02)
    exact_lepton = me*me*ms/(8*np.pi*vh*vh)*(1-4*me*me/ms**2)**1.5
    check('Submuon lepton width agrees with direct two-body expression',
          abs(widths.lepton(ms)/exact_lepton-1) < 3e-15)
    tau_products = [k*k*HBAR/(visible_coefficient(1,.01,widths.unit(ms))*k*k)
                    for k in np.geomspace(1e-11,1e-3,17)]
    check('Portal-only lifetime times kappa squared is invariant',
          np.ptp(tau_products)/np.mean(tau_products) < 1e-14)
    scan = json.loads((ROOT/'SCAN_RESULTS.json').read_text())
    heat = json.loads((ROOT/'HEAT_RESULTS.json').read_text())
    ref = scan['selected'][0]
    recalculated = np.sqrt(HBAR/visible_coefficient(1,.01,widths.unit(ms)))
    check('Original one-second coupling and LZ-screen conflict reproduce',
          abs(recalculated/heat['kappa_for_1s']-1) < 2e-12
          and recalculated/ref['kappa_LZ_shape_matched_ESTIMATE'] > 3e5,
          dict(kappa_one_second=recalculated,
               ratio=recalculated/ref['kappa_LZ_shape_matched_ESTIMATE']))

    profiles_record = json.loads((ROOT/'PROFILE_RESULTS.json').read_text())
    fixed = json.loads((ROOT/'FIXED_YIELD_RESULTS.json').read_text())
    profile_rows = profiles_record['profiles']+[profiles_record['refinement'],
                   profiles_record['old_charge_refinement'], fixed['new_profile'], fixed['refinement']]
    residuals = []
    for p in profile_rows:
        a = np.loadtxt(ROOT/p['profile_file'], delimiter=',', skiprows=1)
        r,f,fp,s,sp = a.T
        integ = lambda u: 4*np.pi*simpson(r*r*u,x=r)
        om, lam = p['Omega'], p['lambda_phi']
        qbar = om*integ(f*f)
        t = integ((fp*fp+sp*sp)/2)
        kinetic = integ(om*om*f*f/2)
        potential = integ(lam/.01*(s*s-1)**2/4+s*s*f*f/2+.1*f**4/4)
        ebar = t+kinetic+potential
        values = dict(Q=p['Q'],lambda_phi=lam,
                      charge=qbar/(.01*p['Q'])-1,
                      energy=ebar/(.1*p['energy_per_v_GeV'])-1,
                      virial=(t+3*(potential-kinetic))/ebar,
                      depletion=(ebar-lam/.01*integ(1-s*s))/ebar)
        residuals.append(values)
    diagnostics['raw_profile_integrals'] = residuals
    check('All archived stationary profiles reproduce charge, energy, virial and depletion from raw arrays',
          max(abs(r[k]) for r in residuals for k in ('charge','energy','virial','depletion')) < 3e-6,
          max(abs(r[k]) for r in residuals for k in ('charge','energy','virial','depletion')))
    refinements = [(next(x for x in profiles_record['profiles'] if x['Q']==1000000 and x['lambda_phi']==.03),profiles_record['refinement']),
                   (next(x for x in profiles_record['profiles'] if x['Q']==39944 and x['lambda_phi']==.021),profiles_record['old_charge_refinement']),
                   (fixed['new_profile'],fixed['refinement'])]
    deltas=[abs(a['energy_per_v_GeV']/b['energy_per_v_GeV']-1) for a,b in refinements]
    check('Larger-domain tighter-tolerance profile energies agree',max(deltas)<2e-6,max(deltas))
    radial=json.loads((ROOT/'RADIAL_RESULTS.json').read_text())
    check('Declared finite-grid radial amplitude Hessians have positive lowest eigenvalues',
          min(x['lowest_fixed_charge_radial_amplitude_eigenvalue'] for x in radial['rows']) > 0.)

    populations = []
    for r in scan['selected']+fixed['rows']:
        pop = r.get('joint_screen_population', r.get('witness_population'))
        if pop is None:
            continue
        q,v,lam=r['Q'],r['v_phi_GeV'],r['lambda_phi']
        energy=r['mass_GeV'];ms=np.sqrt(2*lam)*v;mx=.1*v;md=np.sqrt(.0101)*v
        release=mx*q-energy
        number=release/ms
        bx,bd,bvis=[pop[k] for k in ('hidden_X_branch','hidden_D_branch','visible_branch')]
        # Equal positive and negative daughter charges: their net charge is zero.
        nplus=nminus=bx*number
        nd=2*bd*number
        mass=energy+mx*(nplus+nminus)+md*nd
        yq=pop['required_charge_yield']
        populations.append(dict(Q=q,lambda_phi=lam,v_phi_GeV=v,
                               branch_sum=bx+bd+bvis-1,
                               emitted_energy=number*ms/release-1,
                               late_mass=mass/pop['late_mass_per_object_GeV']-1,
                               net_charge=(q+nplus-nminus)/q-1,
                               matter=yq/q*mass/BDM-1,
                               non_rest_release=release-(mass-energy)))
    diagnostics['population_reconstruction']=populations
    check('Decay branches, charge, source energy and late matter normalization conserve',
          max(abs(x[k]) for x in populations for k in ('branch_sum','emitted_energy','late_mass','net_charge','matter'))<2e-13
          and min(x['non_rest_release'] for x in populations)>-1e-10)
    check('Both new witnesses retain original primordial charge yield',
          max(abs(x['joint_screen_population']['required_charge_yield']/YQ_OLD-1) for x in fixed['rows'])<2e-12
          and all(x['prompt_and_median_projection_screen_overlap'] for x in fixed['rows']))

    detector_checks=[]
    for p,v in [(profile_rows[12],1.),(fixed['new_profile'],fixed['rows'][1]['v_phi_GeV'])]:
        obj=Object(p)
        pop=population(p['Q'],obj.mass(v),v,p['lambda_phi'],1e-9,widths.unit(np.sqrt(2*p['lambda_phi'])*v))
        frac=pop['object_halo_fraction']
        d8,t8,a8=detector(obj,v,frac,points=1601)
        d16,_,_=detector(obj,v,frac,points=6401)
        t16=20*base.bin_integrals(lambda e:frac*obj.rate(e,v),PUB['future_response_xedges'],16)
        a16=t16@PUB['future_response_values']
        errors=dict(Q=p['Q'],lambda_phi=p['lambda_phi'],
                    count_kappa_relative=d8['kappa_LZ_same_mass_count_ESTIMATE']/d16['kappa_LZ_same_mass_count_ESTIMATE']-1,
                    accepted_count_relative=sum(a8)/sum(a16)-1,
                    normalized_zero_momentum_form_factor=float(obj.form(0)-1))
        detector_checks.append(errors)
    diagnostics['detector_refinements']=detector_checks
    check('Recoil quadrature and efficiency grid refine on original and new profiles',
          max(abs(x[k]) for x in detector_checks for k in ('count_kappa_relative','accepted_count_relative','normalized_zero_momentum_form_factor'))<2e-3)

    if cache is not None:
        ys,toys,digest=load_toys(cache)
        controls=[]
        for r,p in zip(fixed['rows'],[next(x for x in profiles_record['profiles'] if x['Q']==39944 and x['lambda_phi']==.021),fixed['new_profile']]):
            obj=Object(p)
            _,_,signals=detector(obj,r['v_phi_GeV'],r['object_halo_fraction'])
            strength,curves=original.curve_grid(signals,ys,toys,2400)
            ul=original.limits(strength,curves,sum(signals))
            pct=np.quantile(1e-9*np.sqrt(ul),[.1586552539,.5,.8413447461])
            err=max(abs(pct/np.array(r['kappa_projected90_percentiles'])-1))
            controls.append(dict(Q=r['Q'],percentile_relative_grid_change=err,
                                 refined_percentiles=pct.tolist(),
                                 upper_screen_still_above_median=bool(pct[1]<r['approximate_LZ_screen_upper'])))
        diagnostics['public_toy_strength_grid_refinement']=controls
        check('Original 1000 public toys refine and retain both approximate overlaps',
              max(x['percentile_relative_grid_change'] for x in controls)<.003
              and all(x['upper_screen_still_above_median'] for x in controls))
    result=dict(scope='Same-author analytic/numerical checks; no independent replication or empirical validation.',
                total_checks=len(checks),passed=sum(x['passed'] for x in checks),
                checks=checks,diagnostics=diagnostics)
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for c in checks:print('PASS' if c['passed'] else 'FAIL',c['name'])
    assert all(x['passed'] for x in checks), 'A deciding control failed; inspect VERIFICATION_RESULTS.json.'


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--cache')
    run(parser.parse_args().cache)

#!/usr/bin/env python3
"""Independent 1D invariant-phase-space integration of the vacuum Born pairs.

Uses the original Codex ground profile. Does not import or execute Claude's
two-momentum/angular quadrature or its profile solver. No medium certification.
"""
import json,sys
from math import comb
from pathlib import Path
import numpy as np
from scipy.integrate import simpson

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
sys.path.insert(0,str(BASE/'grut_core_conversion_20261023/code'))
from scattering import Background,MX,V,G
HBAR=6.582119569e-25
MRED=2.435e18
QREF=200000.
AGE=13.8e9*365.25*86400


def transform(r,source,q):
    out=[]
    for i in range(0,len(q),48):
        out.extend(4*np.pi*simpson((r*r*source)[None,:]*np.sinc(np.outer(q[i:i+48],r)/np.pi),x=r,axis=1))
    return np.array(out)


def free_pairs(bg,d,nr,nq):
    W=(d-2)*bg.omega
    if W<=2:return dict(power=d,open=False,pair_fractional_core_rate_per_second=0.)
    # q=qmax sin(theta): removes the square-root endpoint in invariant phase space.
    r=np.linspace(0,70,nr);F,S=bg.fields(r)
    theta,w=np.polynomial.legendre.leggauss(nq)
    theta=(theta+1)*np.pi/4;w=w*np.pi/4
    qmax=np.sqrt(W*W-4);q=qmax*np.sin(theta)
    ft=transform(r,F**(d-2),q)
    beta=np.sqrt(np.maximum(1-4/(W*W-q*q),0.))
    integrand=q*q*ft*ft*beta*qmax*np.cos(theta)
    integral=np.sum(w*integrand)
    lam=1/MRED**(d-4)
    pref=comb(d,2)**2*lam**2*(V/np.sqrt(2))**(2*d-4)/(8*np.pi**3*MX**3)
    events=pref*integral/HBAR
    rate=(d-2)*events/QREF
    volume=4*np.pi*simpson(r*r*F**(2*d-4),x=r)
    lda=comb(d,2)**2*lam**2*(V/np.sqrt(2))**(2*d-4)*volume*np.sqrt(1-4/W**2)/(4*np.pi*MX**3*HBAR)
    return dict(power=d,open=True,radial_points=nr,phase_points=nq,pair_energy_over_mX=W,
                invariant_q_max=qmax,dimensionless_integral=float(integral),events_per_second=float(events),
                pair_fractional_core_rate_per_second=float(rate),lda_vacuum_fractional_core_rate_per_second=float((d-2)*lda/QREF),
                core_charge_lost_per_event=d-2,escaping_charge_per_event=-2,total_X_violation_per_event=d,
                sampled_q=q[::max(1,nq//96)].tolist(),sampled_integrand=(q*q*ft*ft*beta)[::max(1,nq//96)].tolist())


def main():
    bg=Background();rows=[free_pairs(bg,d,nr,nq) for d in [5,6] for nr,nq in [(6001,96),(12001,192),(16001,384)]]
    ref=json.loads((ROOT/'inputs/CLAUDE_PAIR_REFERENCE.json').read_text())['results']
    old=json.loads((BASE/'grut_cstar_audit_20261024/CHARGE_LEAKAGE_RESULTS.json').read_text())['rows']
    match=json.loads((BASE/'grut_z6_gauge_bridge_20261025/MATCHING_RESULTS.json').read_text())
    summaries=[]
    for d in [5,6]:
        r=[z for z in rows if z['power']==d];fine=r[-1]
        coherent=next(z['initial_fractional_core_charge_loss_rate_per_second_for_c1'] for z in old if z['power']==d and z['radius']==60. and z['potential_scale']==1.)
        pair=fine['pair_fractional_core_rate_per_second'];total=coherent+pair
        summaries.append(dict(power=d,coherent_rate_per_second=coherent,free_Born_pair_rate_per_second=pair,
                              included_channels_rate_per_second=total,pair_to_coherent=pair/coherent,
                              rate_times_age=total*AGE,coefficient_max_for_one_percent=float(np.sqrt(.01/(total*AGE))),
                              S_min_if_abs_c_equals_exp_minusS=float(-np.log(np.sqrt(.01/(total*AGE)))),
                              quadrature_relative_spread=float(np.ptp([z['pair_fractional_core_rate_per_second'] for z in r])/pair),
                              relative_difference_from_review_Born=float(pair/ref[str(d)]['frac_core_rate_free_born']-1)))
    # An analytic Gaussian control checks the radial Fourier normalization.
    r=np.linspace(0,70,14001);q=np.linspace(0,2,101);L=3.
    calc=transform(r,np.exp(-r*r/(2*L*L)),q)
    exact=(2*np.pi)**1.5*L**3*np.exp(-L*L*q*q/2)
    gaussian=np.max(np.abs(calc-exact))/exact[0]
    eff=match['matching']['effective_c6_at_reduced_Planck']
    ground=summaries[1]['included_channels_rate_per_second']*eff*eff
    out=dict(scope='Coherent response plus an independently integrated vacuum Born pair estimate. Not a complete in-medium quantum rate or population lifetime.',
             parent_profile=str(bg.path.relative_to(BASE)),parent_profile_sha256=bg.hash,rows=rows,summaries=summaries,
             controls=dict(Gaussian_Fourier_max_error_over_peak=float(gaussian)),
             reference_age_seconds=AGE,matched_gauge_point=dict(effective_c6=eff,included_channels_rate_per_second=ground,
                 included_channels_rate_times_age=ground*AGE,passes_one_percent_reference_target=ground*AGE<.01),
             missing='Distorted-wave/Bogoliubov pair amplitudes, other multiparticle channels, evolving profiles and population average.')
    (ROOT/'PAIR_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for z in summaries:print('PAIR',z['power'],z['free_Born_pair_rate_per_second'],'BOUND',z['coefficient_max_for_one_percent'],'REVIEW_DIFF',z['relative_difference_from_review_Born'])
    print('GAUGE_POINT_INCLUDED_RATE',ground)


if __name__=='__main__':main()

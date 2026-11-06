#!/usr/bin/env python3
"""Evaluate the fixed H55 dark relic on a declared conditional history."""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid

from production import ARCHIVE, TRACKS, read_columns, sha

S0 = 2891.2
CRIT_OVER_H2 = 1.87834e-29 / 1.7826619216e-24
OMEGA_TARGET = 0.120
E_HIDDEN = 1e-6


def match_release_clock(module, work):
    """Use one interpolant for ln a, H and the derivative of H."""
    source=module.SRC.read_text()
    old='''    herm(TT[i], TT[i + 1], TLA[i], TLA[i + 1], TH[i], TH[i + 1], t, &o->la, NULL);
    herm(TT[i], TT[i + 1], TH[i], TH[i + 1], SH[i], SH[i + 1], t, &o->H, &o->Hd);'''
    new='''    herm(TT[i], TT[i + 1], TLA[i], TLA[i + 1], TH[i], TH[i + 1], t, &o->la, &o->H);
    double d=TT[i+1]-TT[i], s=(t-TT[i])/d;
    o->Hd=((12*s-6)*(TLA[i]-TLA[i+1])+(6*s-4)*d*TH[i]+(6*s-2)*d*TH[i+1])/(d*d);'''
    if source.count(old)!=1:
        raise ValueError('Inherited release interpolation anchor changed')
    module.SRC=work/'release_clock.c'
    module.SRC.write_text(source.replace(old,new))


def vector_readout(assemble, momenta, modes, e, background):
    row=assemble(momenta,modes,e,background)
    scale=momenta**3/(2*np.pi**2)
    first=scale*np.array([r['n1'] for r in modes])
    last=scale*np.array([r['n2'] for r in modes])
    log_q=np.log(momenta)
    row['integrated_readout_change']=float(abs(np.trapezoid(last-first,log_q))/row['N_A'])
    row['absolute_weighted_readout_change']=float(np.trapezoid(abs(last-first),log_q)/row['N_A'])
    return row


def cycle_mean(t,y,width=.025):
    area=cumulative_trapezoid(y,t,initial=0)
    lower=np.maximum(t-width/2,t[0])
    upper=np.minimum(t+width/2,t[-1])
    return (np.interp(upper,t,area)-np.interp(lower,t,area))/(upper-lower)


def table(source, work, kinetic, card, common):
    raw=read_columns(source/'history.txt')
    late=common.read_rec(kinetic/'native')
    ledger=read_columns(kinetic/'native_ledger.txt')
    gap=np.diff(ledger['t'])
    if np.any(gap<0):
        raise ValueError('Kinetic records out of order')
    # The inherited engine can print the final step twice.
    keep=np.r_[gap>0,True]
    late={k:v[keep] for k,v in late.items()}
    ledger={k:v[keep] for k,v in ledger.items()}
    P=card['P']
    q=raw['t']
    tq=np.r_[0.,q[q>0]]
    if abs(tq[-1]-ledger['t'][0])>1e-12:
        raise ValueError('Mismatched coherent and conditional histories')
    def sample(k):
        return np.interp(tq,q,raw[k])
    a_slice=float(np.interp(0.,q,raw['a']))
    ar=sample('a')
    a=np.r_[ar,ledger['a'][1:]]/a_slice
    t=np.r_[tq,ledger['t'][1:]]
    rho=ledger['rho']
    H=np.r_[sample('H'),np.sqrt(rho[1:]/(3*P['mP']**2))]
    R_early=cycle_mean(q,(raw['rho']-3*raw['P'])/P['mP']**2)
    R_late=(rho-3*ledger['pressure'])/P['mP']**2
    R=np.r_[np.interp(tq,q,R_early),R_late[1:]]
    ms2_early=cycle_mean(q,P['g_s']**2*raw['chi2'])
    ms2_late=P['g_s']**2*(P['x0']**2+rho*late['f_chi']/P['m_chi']**2)
    ms2=np.r_[np.interp(tq,q,ms2_early),ms2_late[1:]]
    f_s=np.r_[sample('rho_s')/sample('rho'),(ledger['rho_s']/rho)[1:]]
    f_rad=np.r_[(sample('rho_sm')+sample('rho_dr'))/sample('rho'),
                (late['f_N']+(ledger['rho_sm']+ledger['rho_dr'])/rho)[1:]]
    c=np.r_[sample('s')/P['s_star'],np.sqrt(2*ledger['rho_s'][1:]/ms2_late[1:])/P['s_star']]
    # For the Hermite scale interpolant, H is a quadratic in the
    # fractional step. Test its interior minimum as well as endpoints.
    step=np.diff(t)
    secant=np.diff(np.log(a))/step
    quadratic=3*(H[:-1]+H[1:]-2*secant)
    linear=6*secant-4*H[:-1]-2*H[1:]
    vertex=np.divide(-linear,2*quadratic,out=np.zeros_like(linear),where=quadratic!=0)
    inside=(vertex>0)&(vertex<1)&(quadratic>0)
    lowest=np.min(np.r_[H,(H[:-1]+linear*vertex+quadratic*vertex**2)[inside]])
    if lowest<=0:
        raise ValueError('Nonpositive expansion rate in release interpolant')
    data=np.column_stack((t,np.log(a),H,R,c,np.gradient(c,t),ms2))
    path=work/'history.txt'
    np.savetxt(path,data,fmt='%.17g',header=f'{len(t)} 7',comments='')
    return path,dict(t=t.tolist(),a=a.tolist(),H=H.tolist(),R=R.tolist(),
                     f_s=f_s.tolist(),f_rad=f_rad.tolist()),a_slice,float(max(ms2)**.5)


def run(work, source, kinetic):
    work.mkdir(parents=True,exist_ok=True)
    os.environ['EK_SCRATCH']=str(kinetic)
    os.environ['E2C_SCRATCH']=str(work/'release_work')
    sys.path.insert(0,str(TRACKS/'e2e_kinetic/code'))
    import ek_common as common
    sys.path.insert(0,str(TRACKS/'e2e_chain/code'))
    import e2c_common as E
    import chain
    import relic

    match_release_clock(E.mc,work)
    H55=json.loads(ARCHIVE.read_text())['run']
    card=json.loads((TRACKS/'e2e_card/CANDIDATE_M2.json').read_text())['M2-X']
    c0=card['c0']['value']
    delta=card['holder_phase_delta']['value']
    if delta!=0:
        raise ValueError('Expected the M2-X zero-offset holder')
    Hi=H55['H_GeV']
    sstar=H55['P']['s_star']*Hi
    Gamma_s=1.3e-6*Hi/common.H_REF
    kin=json.loads((kinetic/'summary.json').read_text())
    tab,h,a_slice,bcap=table(source,work,kinetic,H55,common)
    bench=work/'bench.json'
    E.dump(bench,dict(params=dict(H_i_GeV=Hi,s_star_MP=sstar/E.MP_GEV),
                      spectator=dict(Gamma_s_assumed=Gamma_s),history=h))
    original=E.rm.spec_env

    def env(*args,**kwargs):
        setup,gamma=original(*args,**kwargs)
        setup['MFI_SM2TAB']='1'
        setup['MFI_BR']=format((Hi/sstar)**2,'.17g')
        return setup,gamma

    E.rm.spec_env=env
    K=kin['K_total']/a_slice**2
    tag='x11_h55'
    # Match the same holder product at the archived spectator amplitude.
    # Varying the initial spectator value would change the H55 source.
    holder=E.mc.U_PREP
    E.mc.U_PREP=E.mc.MU+c0*(holder-E.mc.MU)
    try:
        release=E.rm.release(tag,K,3.,0.,gs=Gamma_s,sstar=sstar,c0=1.,
                             bg='table:'+str(tab),Tl=1e4,bcap=max(4.,bcap))
    finally:
        E.rm.spec_env=original
        E.mc.U_PREP=holder
    assemble=E.vm.assemble
    E.vm.assemble=lambda *args: vector_readout(assemble,*args)
    try:
        E.vm.main(tag)
    finally:
        E.vm.assemble=assemble
    vectors=E.load(E.ENG/'vectors'/f'vectors_{tag}.json')
    vector=next(r for r in vectors['by_e'] if r['e']==E_HIDDEN)
    E.dump(E.ENG/'chain'/f'chain_{tag}.json',E.cm.analyse(tag))
    gas=chain.full(tag,bench,c0,delta)['by_e']['cL0_vanish_at_MY']['e1e-6']
    charge=abs(release['result']['N_X'])
    entropy=kin['visible_entropy_comoving_Hi3']/a_slice**3
    YQ=charge/entropy
    # The population model receives the calculated YQ. Its abundance
    # calibration is deliberately bypassed.
    population=relic.relic(gas['S_per_Q'],H_i=Hi,e=E_HIDDEN,Y_Q=YQ)
    M=population['M_central']
    Yreq=OMEGA_TARGET*CRIT_OVER_H2/(relic.MX*S0)
    out=dict(scope='conditional relic from dephased production and inherited hidden release and population models',
             H55_scale_GeV=Hi,holder_c0=c0,holder_offset_degrees=delta,e_hidden=E_HIDDEN,
             holder_strength_implementation='scaled coupling at fixed archived spectator amplitude',
             source_time=kin['source_time'],K_total=K,K_visible=kin['K_visible']/a_slice**2,
             s_comoving_visible_Hi3=entropy,release=release['result'],
             ecrit_used_for_certification=False,hidden_S_per_Q=gas['S_per_Q'],
             hidden_heat_terms={k:gas[k] for k in ('R_total','R_vec','R_Y','R_casc','R_leak','dR_grav')},
             vector_readout={k:vector[k] for k in ('max_readout_change','integrated_readout_change',
                                                  'absolute_weighted_readout_change')},
             X_mass_GeV=relic.MX,Y_charge=YQ,population=population,
             total_relic_mass_per_visible_entropy_GeV=float(M*relic.MX*YQ),
             conditional_Omega_h2=float(M*relic.MX*YQ*S0/CRIT_OVER_H2),
             target_Omega_h2=OMEGA_TARGET,target_Y_if_all_free_X=float(Yreq),
             source_handoff_passed=kin['initial_splice']['passed'],
             release_clock_source_sha256=sha(E.mc.SRC),
             source_table_sha256=sha(tab),accepted_abundance=False)
    (work/'summary.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(f'{work.name}: conditional Omega h^2 {out["conditional_Omega_h2"]:.6g}, '
          f'YQ={YQ:.6g}, S/Q={gas["S_per_Q"]:.6g}',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('source',type=Path)
    parser.add_argument('kinetic',type=Path)
    args=parser.parse_args()
    run(args.work.resolve(),args.source.resolve(),args.kinetic.resolve())


if __name__=='__main__':
    main()

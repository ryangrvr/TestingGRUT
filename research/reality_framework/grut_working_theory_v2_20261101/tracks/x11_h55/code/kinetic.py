#!/usr/bin/env python3
"""Conditional kinetic continuation of the native H55 occupation spectra."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from production import TRACKS, ARCHIVE, read_columns, sha

sys.path.insert(0,str(TRACKS/'x4_closure/code'))
from kinetic_engine import source as x4_source
from kinetic_replay import input_tables


def compile_engine(work, step_factor, legacy_expansion):
    source = x4_source(step_factor)
    old = '%.17g %.17g %.17g %.17g\\n",\n            t, a, rho_s, rho_dr, rho_h + rho_th, rho, exchange_error,\n            budget_error, occ_min, occ_max);'
    new = '%.17g %.17g %.17g %.17g %.17g\\n",\n            t, a, rho_s, rho_dr, rho_h + rho_th, rho, exchange_error,\n            budget_error, occ_min, occ_max, pressure());'
    header = 'budget_error f_min f_max\\n'
    if source.count(old)!=1 or source.count(header)!=1:
        raise ValueError('Inherited pressure ledger anchor changed')
    source=source.replace(old,new).replace(header,'budget_error f_min f_max pressure\\n')
    if not legacy_expansion:
        clock='double aold = a, anew = noexp ? a : a * exp(H * dt);'
        if source.count(clock)!=1:
            raise ValueError('Inherited expansion anchor changed')
        source=source.replace(clock,'''double H0 = sqrt(exchanged / (3*mP*mP));
        double Hdot = -(exchanged + pressure_before) / (2*mP*mP);
        double aold = a, anew = noexp ? a : a * exp(dt * (H0 + 0.5*dt*Hdot));''')
    src = work/'engine.c'
    src.write_text(source)
    binary = work/'engine'
    subprocess.run(['gcc','-std=gnu11','-O3','-march=native',str(src),'-lm','-o',str(binary)],check=True)
    return binary


def native_state(history, spectra, card, common, added_bins):
    last = {k:float(v[-1]) for k,v in history.items()}
    a = last['a']
    masses = np.array([ch['M'] for ch in card['channels']])
    q_source = spectra['q'][(spectra['species']==1)&(spectra['h']==-1)]
    w_source = spectra['w'][(spectra['species']==1)&(spectra['h']==-1)]
    dq = 2*q_source[0]
    q_end = q_source[-1]+dq/2
    extra, weights = common.grid_log(q_end,3e7,added_bins)
    q, w = np.r_[q_source,extra], np.r_[w_source,weights]
    states, moments = [], []
    for i,M in enumerate(masses,1):
        selected = spectra['species']==i
        native = spectra['f_free'][selected].reshape(2,-1)
        expected_q = np.tile(q_source,2)
        if not np.array_equal(spectra['q'][selected],expected_q):
            raise ValueError('Unexpected native helicity or momentum layout')
        if native.min() < -1e-10 or native.max() > 1+1e-10:
            raise ValueError('Source occupations outside Pauli interval')
        occ = np.pad(native,((0,0),(0,len(extra))))
        states.append(occ.ravel())
        f = occ.sum(axis=0)
        E = np.hypot(q/a,M)
        moments.append([float(w@f/a**3),float((w*E)@f/a**3),float((w*q*q/(3*E*a*a))@f/a**3)])
    m = np.array(moments)
    expected = np.array([[last[f'n_{i}'],last[f'rho_{i}'],last[f'p_free_{i}']] for i in range(1,5)])
    fidelity = float(np.max(abs(m/expected-1)))
    if fidelity>1e-10:
        raise ValueError('Native energy, pressure or number projection failed')
    field_energy = last['rho']-m[:,1].sum()-last['rho_s']-last['rho_sm']-last['rho_dr']
    # A declared harmonic partition for the two parent reservoirs. This
    # scales their positive harmonic energies to the evolved total and
    # retains the full native fermion state projected onto its occupations.
    P = card['P']
    harmonic = []
    for j,name in enumerate(('sigma','chi')):
        phase_speed = last[f'drive_{j+1}']*2*masses[j]*card['channels'][j]['f']
        delta2 = (last[name]-P[('s0','x0')[j]])**2+last[f'rms_{name}']**2
        harmonic.append((phase_speed**2+P[('m_sig','m_chi')[j]]**2*delta2)/2)
    parents = field_energy*np.array(harmonic)/sum(harmonic)
    if min(parents)<0:
        raise ValueError('Negative harmonic-parent energy')
    state = dict(t=last['t'],a=a,mP=P['mP'],rho_sig=parents[0],rho_chi=parents[1],t_end=1e8,
                 rho_h=0.,n_h=0.,rho_th=last['rho_sm'])
    after_p = (last['rho_sm']+last['rho_dr'])/3+m[:,2].sum()
    sampled_work = sum(history[f'work_{i}'] for i in range(1,5))/history['rho']
    tail = history['t']>=last['t']-.02
    splice = dict(time=last['t'],parent_energy_partition=parents.tolist(),
                  native_moment_max_relative_error=fidelity,
                  fermion_pressure_jump_over_density=float((m[:,2].sum()-sum(last[f'p_{i}'] for i in range(1,5)))/last['rho']),
                  complete_pressure_jump_over_density=float((after_p-last['P'])/last['rho']),
                  sampled_residual_work_over_density=float(sampled_work[-1]),
                  max_abs_sampled_residual_work_last_window=float(max(abs(sampled_work[tail]))))
    return state,q,w,states,splice


def run(work, source_path, extra_bins, step_factor, later_n2, legacy_expansion):
    work.mkdir(parents=True,exist_ok=True)
    os.environ['EK_SCRATCH']=str(work)
    sys.path.insert(0,str(TRACKS/'e2e_kinetic/code'))
    import ek_common as common

    card=json.loads(ARCHIVE.read_text())['run']
    history=read_columns(source_path/'history.txt')
    spectra=read_columns(source_path/'spectra.txt')
    st,q,w,occ,splice=native_state(history,spectra,card,common,extra_bins)
    source={k:float(v[-1]) for k,v in history.items()}
    if later_n2:
        seesaw=json.loads((TRACKS/'e2e_sm/results/seesaw_H1.3e11.json').read_text())
        row=next(row for row in seesaw['scenarios']['S1_mt050']['rows'] if row['delta']==0)
        tilde=row['mtilde_eV']
    else:
        tilde=[.05]*4
    spec,mpar=common.species('HFa',card['H_GeV'],[x*1e-9 for x in tilde])
    gamma_s=1.3e-6*card['H_GeV']/common.H_REF
    inp=work/'native.in'
    common.write_input(inp,st,spec,mpar,common.sm_coeffs(),dict(coherent=0),[(q,w)]*4,occ,
                       tables=input_tables(),dshell=-1e5)
    with inp.open('a') as f:
        f.write(f"{source['rho_s']:.17g} {source['rho_dr']:.17g} {gamma_s:.17g}\n")
    binary=compile_engine(work,step_factor,legacy_expansion)
    prefix=work/'native'
    process=subprocess.run([str(binary),'run',str(inp),str(prefix)],text=True,
                           capture_output=True,timeout=2400)
    if process.returncode:
        raise RuntimeError(process.stderr)
    rec=common.read_rec(prefix)
    ledger=read_columns(work/'native_ledger.txt')
    if not np.allclose(ledger['t'],rec['t'],rtol=6e-11,atol=0):
        raise ValueError('Kinetic diagnostic times differ')
    a,t=ledger['a'],ledger['t']
    rho=ledger['rho']
    K_total=a*a*np.sqrt(rho/(3*st['mP']**2))
    K_visible=a*a*np.sqrt(ledger['rho_sm']/(3*st['mP']**2))
    eq=rec['f_th']>0.5
    i_rh=int(np.argmax(eq))
    if not np.any(eq):
        raise ValueError('Thermal SM fails to dominate in conditional history')
    temperature_comoving=(30*ledger['rho_sm'][-1]*a[-1]**4/(np.pi**2*common.GSTAR))**0.25
    entropy=(2*np.pi**2/45)*common.GSTAR*temperature_comoving**3
    splice['energy_gap']=float(rho[0]/source['rho']-1)
    splice['Hubble_gap']=float(K_total[0]/a[0]**2/source['H']-1)
    splice['passed']=bool(max(abs(splice[k]) for k in ('energy_gap','complete_pressure_jump_over_density',
                                                       'max_abs_sampled_residual_work_last_window'))<1e-6)
    idx=sorted({0,len(a)-1,*(int(np.argmin(abs(t-x))) for x in (1,3,10,100,1e3,1e4,1e5,1e6,1e7))})
    result=dict(scope='conditional dephased kinetic continuation from native occupations',
                source_time=st['t'],extra_momentum_bins=extra_bins,native_momentum_bins=len(q)-extra_bins,
                mtilde_eV=tilde,step_factor=step_factor,legacy_expansion=legacy_expansion,initial_splice=splice,
                late_time=float(t[-1]),K_total=float(K_total[-1]),K_visible=float(K_visible[-1]),
                visible_entropy_comoving_Hi3=float(entropy),
                rho_DR_over_visible=float(ledger['rho_dr'][-1]/ledger['rho_sm'][-1]),
                Delta_Neff_if_standard_entropy=float((43/7)*(10.75/common.GSTAR)**(1/3)*ledger['rho_dr'][-1]/ledger['rho_sm'][-1]),
                T_SM_half_GeV=float(rec['T_th'][i_rh]*card['H_GeV']),
                T_SM_half_time=float(t[i_rh]),
                K_visible_late_drift=float(K_visible[-1]/np.interp(t[-1]/10,t,K_visible)-1),
                max_kinetic_budget_error=float(max(ledger['budget_error'])),
                max_kinetic_exchange_error=float(max(ledger['exchange_error'])),
                occupation_range=[float(min(ledger['f_min'])),float(max(ledger['f_max']))],
                source_history_sha256=sha(source_path/'history.txt'),
                source_spectra_sha256=sha(source_path/'spectra.txt'),
                kinetic_input_sha256=sha(inp),
                samples=[dict(t=float(t[i]),a=float(a[i]),H=float(K_total[i]/a[i]**2),rho=float(rho[i]),
                              P=float(ledger['pressure'][i]),K_visible=float(K_visible[i]),
                              f_thermal=float(rec['f_th'][i])) for i in idx],
                accepted_abundance=False)
    (work/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f"{work.name}: K_visible={K_visible[-1]:.9g}, "
          f"energy gap={splice['energy_gap']:.3g}, pressure gap={splice['complete_pressure_jump_over_density']:.3g}",
          flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('source',type=Path)
    parser.add_argument('--extra-bins',type=int,default=2400)
    parser.add_argument('--step-factor',type=float,default=1.)
    parser.add_argument('--later-n2',action='store_true')
    parser.add_argument('--legacy-expansion',action='store_true')
    args=parser.parse_args()
    run(args.work,args.source,args.extra_bins,args.step_factor,args.later_n2,args.legacy_expansion)


if __name__=='__main__':
    main()

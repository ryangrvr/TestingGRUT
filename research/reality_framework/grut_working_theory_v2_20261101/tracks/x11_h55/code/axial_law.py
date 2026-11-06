#!/usr/bin/env python3
"""Apply the H55 axial CP matching postulate to a fixed seesaw texture."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

TRACK = Path(__file__).resolve().parents[1]
TRACKS = TRACK.parent
sys.path.insert(0, str(TRACKS/'e2e_sm/code'))
from sm_seesaw import V_EW, fS, fV, mnu_loop

SEESAW = TRACKS/'e2e_sm/results/seesaw_H1.3e11.json'
H55 = TRACKS/'x2_strip/results/lower_H_lattice.json'
ETA_OBS = 6.12e-10


def evaluate(real_yukawa, M, phase, parity):
    roots = np.where(parity > 0, 1.+0j, 1j)
    Y = real_yukawa * (roots*np.exp(-0.5j*phase))
    grams = Y.conj().T @ Y
    widths = np.real(np.diag(grams))*M/(8*math.pi)
    I = np.imag(grams**2)
    norm = abs(grams)**2
    normed = np.divide(I,norm,out=np.zeros_like(I),where=norm>1e-30)
    rule = parity[:,None]*parity[None,:]*np.sin(phase[:,None]-phase[None,:])
    invariant_error = float(np.max(abs(normed-rule)[norm>1e-30]))
    epsilon = np.zeros(len(M))
    for i in range(len(M)):
        for j in range(len(M)):
            if i==j:
                continue
            x = (M[j]/M[i])**2
            epsilon[i]-=I[i,j]*(fV(x)+fS(x))/(8*math.pi*grams[i,i].real)
    mD = V_EW*Y
    mnu = -(mD/M) @ mD.T
    masses = np.sort(np.linalg.svd(mnu,compute_uv=False))*1e9
    tilde = np.sum(real_yukawa**2,axis=0)*V_EW**2/M*1e9
    trace = abs(np.trace(mnu))*1e9
    trace_formula = abs(np.sum(tilde*parity*np.exp(-1j*phase)))
    return dict(parity=parity.tolist(),mtilde_eV=tilde.tolist(),Gamma_rest_GeV=widths.tolist(),
                normalized_CP_invariants=normed.tolist(),epsilon=epsilon.tolist(),
                light_masses_eV=masses.tolist(),sum_light_masses_eV=float(masses.sum()),
                squared_mass_splittings_eV2=(masses[1:]**2-masses[0]**2).tolist(),
                mass_sum_trace_bound_eV=float(trace),trace_formula_error_eV=float(abs(trace-trace_formula)),
                invariant_error=invariant_error,
                repository_loop_estimate=mnu_loop(mD,M))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    source = json.loads(SEESAW.read_text())
    case = json.loads(H55.read_text())['run']
    baseline = next(row for row in source['scenarios']['S1_mt050']['rows'] if row['delta']==0)
    mass = np.array([s['M'] for s in case['channels']])*case['H_GeV']
    phase = np.array([case['P'][('s0','x0')[s['w']]]/s['f'] for s in case['channels']])
    parity = np.array(source['ETA'],int)
    real = (np.array(baseline['R_re'])+np.array(baseline['R_im']))
    real *= np.sqrt(np.array(source['light_masses_eV'])*1e-9)[:,None]*np.sqrt(mass)[None,:]/V_EW
    benchmark_tilde = np.array(baseline['mtilde_eV'])
    widths = np.array([sp['Grest'] for sp in case['channels']])*case['H_GeV']
    tilde_native = 8*math.pi*widths*V_EW**2/mass**2*1e9
    native = real*np.sqrt(tilde_native/benchmark_tilde)
    literal = evaluate(native,mass,phase,np.ones(4,int))
    ext = evaluate(native,mass,phase,parity)
    ext_later = evaluate(real,mass,phase,parity)
    zero = evaluate(real,mass,np.zeros_like(phase),parity)
    common = evaluate(real,mass,phase+0.61,parity)
    conj = evaluate(real,mass,-phase,parity)
    result = dict(
        scope='candidate axial CP matching at tree level with one fixed real texture',
        H55_scale_GeV=case['H_GeV'], sterile_masses_GeV=mass.tolist(),
        vacuum_fields_H55=[case['P']['s0'],case['P']['x0']],
        mass_phase_radians=phase.tolist(),
        reduced_mass_phase_radians=np.angle(np.exp(1j*phase)).tolist(),
        eps_convention='inherited sm_seesaw unflavoured sign convention',
        literal_positive_parity=literal,
        inherited_parity_native_widths=ext,
        inherited_parity_later_widths=ext_later,
        tests=dict(common_phase_invariant_error=float(np.max(abs(
            np.array(common['normalized_CP_invariants'])-ext_later['normalized_CP_invariants']))),
                   common_phase_mass_error_eV=float(np.max(abs(
            np.array(common['light_masses_eV'])-ext_later['light_masses_eV']))),
                   conjugate_decay_epsilon_error=float(np.max(abs(
            np.array(conj['epsilon'])+ext_later['epsilon']))),
                   zero_phase_max_epsilon=float(max(map(abs,zero['epsilon']))),
                   zero_phase_mass_error_eV=float(np.max(abs(
            np.array(zero['light_masses_eV'])-source['light_masses_eV']))),
                   native_width_reproduction_error=float(np.max(abs(np.array(literal['Gamma_rest_GeV'])/widths-1)))),
        Planck_BAO_mass_sum_limit_eV=0.12, benchmark_eta_B=ETA_OBS,
        input_sha256={str(p.relative_to(TRACKS)):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (SEESAW,H55)},
        accepted_GRUT_law=False)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for key in ('literal_positive_parity','inherited_parity_native_widths','inherited_parity_later_widths'):
        row=result[key]
        print(key,f"mass sum {row['sum_light_masses_eV']:.6g} eV, trace bound "
              f"{row['mass_sum_trace_bound_eV']:.6g} eV, epsilon {row['epsilon']}",flush=True)


if __name__=='__main__':
    main()

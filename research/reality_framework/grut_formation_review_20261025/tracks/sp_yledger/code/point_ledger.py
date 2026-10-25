#!/usr/bin/env python3
"""sp_yledger (Claude, own code): the leading 373d37e point and its Y/Sigma/vector debris ledger.

1. Independent release (cosmic time, physical fields, DP5(4); Codex: conformal fields, DOP853) of
   H_i = 1.1211167304108865e13 GeV, c_H = 1, delta = pi/2, to a = 3000; tolerance and preparation controls.
2. Late Y excitation: (i) Codex-style single-mass number J_Y = E_exc/M; (ii) our projection onto the two real
   Majorana eigenmodes, N_Y = sum_j E_j/omega_j, with heavy/light split; Sigma radial number N_sig.
3. Decay ledger [D/E]: widths of Y_L -> 3X, Y_H -> Y_L + A_L (Goldstone equivalence), sigma -> A_L A_L,
   A -> X Xbar, Y -> XX on the condensate; decay epochs; pairs per net charge; deposited energies using the
   exponential-decay average of a(t) in radiation; no energy is counted twice (E_X, E_Y,exc, E_sig are
   disjoint pieces of the Hamiltonian; the A energy is taken out of the Y and sigma budgets).
Usage: python3 point_ledger.py <out.json>     (~1 min)
"""
import json, sys
from multiprocessing import Pool
import numpy as np
from ymodel import initial, run, entropy_inst, MRED, GSTAR, YQ, MY, HV, MU, KAP, LS, E_GAUGE, MXL, BDM

H = 1.1211167304108865e13
CODEX = dict(eps=-0.24066124082075235, NX=-238.64867701011931, EX=741.6498820668295, JY=13.269029096898693,
             a_dec=172230.3091474945, EY_over_EX=274.85203825803666, pairs=0.0834010223509602,
             pair_limit=0.0049296413205678555, eps_energy_equiv=0.0035554862510437835, Gamma=377.9476691390566)
TEND = 0.5 * 3000**2


def job(args):
    label, z0, rtol, atol = args
    r = run(H, z0, t_end=TEND, rtol=rtol, atol=atol, nsamp=120)
    return label, r


def main():
    z0, m = initial(H, 1.0, np.pi / 2)
    zs = z0.copy(); zs[2] = zs[2] * (1 + 2.08e-7)          # radial preparation control (|Sigma| x (1+2e-7))
    zp = z0.copy(); zp[2] = zp[2] * np.exp(1e-9j)          # source phase +1e-9 rad
    jobs = [('nominal', z0, 1e-12, 1e-15), ('rtol1e-13', z0, 1e-13, 1e-16), ('rtol1e-11', z0, 1e-11, 1e-14),
            ('Sigma_radius_x(1+2.08e-7)', zs, 1e-12, 1e-15), ('phase+1e-9', zp, 1e-12, 1e-15)]
    with Pool(4) as p:
        res = dict(p.map(job, jobs))
    nom = res['nominal']
    W = nom['windows']
    keys = ['eps', 'NX', 'EX', 'JYcodex', 'NYH', 'NYL', 'Nsig', 'NY', 'Neff', 'BoverA']
    win = {k: {c: W[k][c] for c in keys} for k in W}
    w6, wl = W['600-1000'], W['2000-3000']
    controls = {lab: dict(end=r['end'], w600=dict(eps=r['windows']['600-1000']['eps'],
                                                    JY=r['windows']['600-1000']['JYcodex'],
                                                    NY=r['windows']['600-1000']['NYH'] + r['windows']['600-1000']['NYL'],
                                                    Nsig=r['windows']['600-1000']['Nsig']))
                for lab, r in res.items()}

    # ---- decay ledger (GeV unless stated) ----
    Svac = HV / np.sqrt(2)
    mH = np.sqrt(MY**2 + 2 * MU * Svac); mL = np.sqrt(MY**2 - 2 * MU * Svac)
    msig = np.sqrt(2 * LS) * HV
    mA = E_GAUGE * 6 * HV                                   # = 6 e h
    G3_M = 3 * KAP**2 * MY / (256 * np.pi**3)               # Codex convention
    G3_L = 3 * KAP**2 * mL / (256 * np.pi**3)               # light real eigenstate: XXX and XbXbXb, each half
    pHL = (mH**2 - mL**2) / (2 * mH)
    G_HL = (np.sqrt(2) * MU)**2 * pHL / (8 * np.pi * mH**2)  # Y_H -> Y_L + G (= A_L), coupling sqrt2 mu G Y1 Y2
    G_sig = msig**3 / (32 * np.pi * HV**2)                   # sigma -> G G
    G_A = E_GAUGE**2 * mA / (48 * np.pi) * (1 - 4 * MXL**2 / mA**2)**1.5
    a_eq = lambda G: np.sqrt(H / G)                         # H = Gamma in radiation (Codex convention)
    a_avg = lambda G: np.sqrt(2 * H / G) * 0.886226925452758  # <a> over exp(-Gamma t), a = sqrt(2 H_i t)
    # Y -> XX with one condensate leg: Gamma2/H constant in radiation; X amplitude from E_X (quartic)
    Xconf = (CODEX['EX'] / 1e-3)**0.25
    G2_over_H = 9 * KAP**2 * (Xconf * H)**2 / (8 * np.pi * MY) / H
    NX = wl['NX']; NY = wl['NYH'] + wl['NYL']; NYH_late = wl['NYH']; NYH_early = W['100-300']['NYH']
    Nsig = wl['Nsig']
    pair_limit = 0.01 * BDM / (2 * MXL * YQ)
    f_per_pair = 2 * MXL * YQ / BDM
    pairs_Y_codex = 1.5 * w6['JYcodex'] / abs(w6['NX'])
    pairs_Y = 1.5 * NY / abs(NX)
    pairs_sig = 2 * Nsig / abs(NX)
    pairs_AY = (NYH_late / abs(NX), NYH_early / abs(NX))
    tot = (pairs_Y + pairs_sig + pairs_AY[0], pairs_Y + pairs_sig + pairs_AY[1])
    # deposited conformal energies (H_i^4 a^4 rho units), compare with E_X
    EX = wl['EX']
    EY_codex = w6['JYcodex'] * (MY / H) * a_eq(G3_M) / w6['EX']
    EY_dep = NY * (mL / H) * a_avg(G3_L) / EX
    EA_unit = (mA / H) * a_avg(G_A) / EX                     # per A quantum, non-relativistic at decay
    EA_dep = ((2 * Nsig + NYH_late) * EA_unit, (2 * Nsig + NYH_early) * EA_unit)
    eps = abs(wl['eps'])
    eq = lambda R: eps / (1 + R)**0.75
    out = dict(
        point=dict(H_GeV=H, cH=1.0, delta=np.pi / 2, minimum_over_H=[m['x'], m['y'], m['s']], grad_norm=m['grad_norm']),
        windows=win, controls=controls,
        reproduction_vs_codex_600_1000=dict(eps=(w6['eps'], CODEX['eps']), NX=(w6['NX'], CODEX['NX']),
                                            EX=(w6['EX'], CODEX['EX']), JY=(w6['JYcodex'], CODEX['JY'])),
        masses_GeV=dict(m_YH=mH, m_YL=mL, m_sigma=msig, m_A=mA, Delta_M=mH - mL),
        widths_GeV=dict(Y_3X_codex_M=G3_M, YL_3X=G3_L, YH_to_YL_A=G_HL, sigma_to_AA=G_sig, A_to_XXbar=G_A,
                        Y_to_XX_condensate_over_H=G2_over_H),
        a_Gamma_eq_H=dict(Y_codex=a_eq(G3_M), YL=a_eq(G3_L), YH_to_YL=a_eq(G_HL), sigma=a_eq(G_sig), A=a_eq(G_A)),
        a_mean_decay_radiation=dict(YL=a_avg(G3_L), A=a_avg(G_A)),
        T_vis_at_A_decay_GeV=entropy_inst(H)[1] / a_avg(G_A),
        numbers_late=dict(NX=NX, NY_modes=NY, NY_heavy_late=NYH_late, NY_heavy_100_300=NYH_early,
                          JY_codex_style=wl['JYcodex'], N_sigma=Nsig, NY_charge=wl['NY'],
                          three_NY_over_NX=3 * wl['NY'] / NX),
        pairs_per_net_charge=dict(codex=pairs_Y_codex, Y_daughters=pairs_Y, sigma_vectors=pairs_sig,
                                  YH_vectors_range=pairs_AY, total_range=tot, one_percent_allowance=pair_limit,
                                  ratio_codex=pairs_Y_codex / pair_limit,
                                  ratio_total_range=(tot[0] / pair_limit, tot[1] / pair_limit),
                                  DM_fraction_unthermalized_range=(tot[0] * f_per_pair, tot[1] * f_per_pair),
                                  X_remainder_frozen_f=1.01427 * (1 / eps - 1)),
        energies_over_EX=dict(Y_codex_convention=EY_codex, Y_deposited=EY_dep, A_per_quantum=EA_unit,
                              A_deposited_range=EA_dep,
                              total_injected_range=(EY_dep + EA_dep[0], EY_dep + EA_dep[1])),
        energy_equivalent_eps=dict(codex=CODEX['eps_energy_equiv'], Y_only=eq(EY_dep),
                                   Y_plus_A_range=(eq(EY_dep + EA_dep[0]), eq(EY_dep + EA_dep[1]))),
        labels=dict(numbers='[C] own release', widths='[D] tree level, vacuum masses', ledger='[E] conditional'))
    json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ('windows',)}, indent=1, default=float))


if __name__ == '__main__':
    main()

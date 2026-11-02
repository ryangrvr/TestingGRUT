#!/usr/bin/env python3
"""inf_potential SKEPTIC (Claude): control scan of the M1r release engine (rpinf.c MODE 0 = rpdq.c, the quadratic
oscillator started at rest at a = 1 with amplitude sqrt(6) H M_P/m_I and mu_eff = 10i there) versus m_I/H_i, Gamma_I = 0.5.
Question: is the large homogeneous Y_H remnant of the MODE-4 benchmarks a property of the inflaton FREQUENCY (a
non-adiabatic / Y_H-resonant release once m_I is comparable to the Y masses m_L = 22.25, m_H = 36.12 H_i), independent of
the potential shape?  Output ../results/sk_mi_scan.json.  Labels: [C] engine; interpretation [A]."""
import json, math, sys
import numpy as np
import sk_engine as se
from pathlib import Path
OUT = Path(__file__).resolve().parent.parent / 'results'

def run_b(b, tend=40.0):
    p = se.fc.params(0.5); z = se.fc.stationary(se.fc.U_PREP, p, R=12.0)['z']
    eta2 = se.fc.ETA2; I0 = math.sqrt(6.0) / b
    G = (se.fc.U_PREP - p['mu']) / (I0 / math.sqrt(eta2))
    args = se.common_args(0, b, 0.5, eta2, p, G, tend, 400) + se.spec_args(z)
    import subprocess
    r = subprocess.run([str(se.SP / 'rpinf')] + args, capture_output=True, text=True)
    S = se.parse_S(r.stdout); late = S[S[:, 0] >= 0.6 * tend]
    K = float(np.mean(late[:, 1] ** 2 * late[:, 2]))
    return dict(m_I=b, NX=float(np.mean(late[:, 4])), EX=float(np.mean(late[:, 8])), eps=float(np.mean(late[:, 9])),
                NYH=float(np.mean(late[:, 10])), NYL=float(np.mean(late[:, 11])), Nsig=float(np.mean(late[:, 12])), K=K)


def sq_ledger(r, e=1e-6, yl=True):
    # dq_heat ledger on this release (static-slice exit factor f_eff = 1; vectors N_A ∝ K^3/2, deposits ∝ K^1/2)
    import sk_reledger as L
    K = r['K']; row = L.rows[e]
    NA = (row['N_A_k_le_1'] + row['N_A_k_ge_1']) * (K / L.KM) ** 1.5
    D = NA * L.depA(e, K) + (r['NYL'] + r['NYH']) * (L.depYL(K) if yl else 0.0) + (r['NYH'] + 2 * r['Nsig']) * L.depA(e, K)
    Ec = (abs(r['NX']) / L.C_HART) ** (4 / 3)
    return L.SQ((r['EX'] + D) / Ec)

if __name__ == '__main__':
    se.build()
    rows = []
    for b in (1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 11, 15, 20, 25, 30, 36, 42.37, 50, 64, 90, 128):
        r = run_b(b)
        r['SQ_e1e-6'] = sq_ledger(r); r['SQ_e1e-5'] = sq_ledger(r, 1e-5); r['SQ_e1e-6_noYLdeposit'] = sq_ledger(r, yl=False)
        rows.append(r)
        print('m_I=%6.2f K=%.3f NX=%9.1f EX=%9.4g eps=%.4f NYH=%9.4g NYL=%8.4g Nsig=%.3g | S/Q(1e-6)=%.2f S/Q(1e-5)=%.2f S/Q(1e-6, no Y_L deposit)=%.2f' % (
            r['m_I'], r['K'], r['NX'], r['EX'], r['eps'], r['NYH'], r['NYL'], r['Nsig'], r['SQ_e1e-6'], r['SQ_e1e-5'], r['SQ_e1e-6_noYLdeposit']), flush=True)
    json.dump(dict(method=__doc__, rows=rows), open(OUT / 'sk_mi_scan.json', 'w'), indent=1)

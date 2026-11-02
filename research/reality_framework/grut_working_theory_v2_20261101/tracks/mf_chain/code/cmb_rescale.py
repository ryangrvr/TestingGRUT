#!/usr/bin/env python3
"""mf_chain: CMB side of the H_i iteration [C on the inherited capacity definition, D for r and N_*].

Rescaling the hybrid's V0 at fixed (m_chi, m_sig)/H_i and fixed kappa = mu_R^2/(g 2 sigma_c)^2 = 10 raises the waterfall's
X = g^2 sigma_c = m_sig m_chi^2/(4 sqrt3 M_P) in H_i units proportionally to H_i, so the chi one-loop term weighs more on
the clock.  The e-fold capacity N_max (largest number of e-folds before sigma_c for which the retuned tree clock (c1, m2)
meets A_s and n_s with a monotonic potential, pivot within 0.2 M_P) is computed with the mf_hybrid skeptic's
sk_structure.nmax (imported READ-ONLY; its module constant mP = M_P/H is patched per H), regression at the design H
against the skeptic's (250, 600, kappa = 10) value 74.3.  N_* comes from mf_hybrid's Nstar_from_K (read-only) with the
history's own K; r = 2 H_*^2/(pi^2 A_s M_P^2) with H_* = H_i (1 + 2.35e-6, the bench's exit value).
Usage: python3 -B cmb_rescale.py H_GeV:K:ns [...]   -> results/cmb_rescale.json"""
import json, math, sys
import mfc_common as M

sys.path.insert(0, str(M.SKC))
sys.path.insert(0, str(M.TRK / 'mf_hybrid' / 'code'))
import sk_common as SC                                   # noqa: E402 read-only
import sk_structure as SS                                # noqa: E402 read-only
import mfh_common as MH                                  # noqa: E402 read-only

AS = 2.1e-9


def capacity(H, kappa=10.0, ns=0.9649):
    mP = M.MP_GEV / H
    SC.mP = mP; SS.mP = mP
    best = (0.0, None)
    P = SS.params(250.0, 600.0, kappa)
    import numpy as np
    for dl in np.geomspace(1e-4, 0.2, 140):
        N = SS.nefolds(P, dl, ns)
        if N > best[0]:
            best = (N, dl)
    return dict(H_GeV=H, mP=mP, X=P['g'] ** 2 * P['sc'], Nmax=best[0], dl_at_Nmax=best[1], sigma_c_MP=P['sc'] / mP)


def main(args):
    out = dict(method=__doc__, rows=[])
    reg = capacity(M.H_REF)
    out['regression_design'] = reg
    print('regression at 1.3e11: Nmax %.2f (skeptic 74.3), X %.4f' % (reg['Nmax'], reg['X']))
    for a in args:
        H, K, ns = (float(x) for x in a.split(':'))
        c = capacity(H, 10.0, ns)
        Hs = H * (1 + 2.35e-6)
        Nst = MH.Nstar_from_K(Hs, H, K)
        r = 2 * Hs ** 2 / (math.pi ** 2 * AS * M.MP_GEV ** 2)
        c.update(K=K, ns_target=ns, Nstar=Nst, r=r, capacity_ok=c['Nmax'] >= Nst)
        out['rows'].append(c)
        print('H %.4e K %.3f ns %.4f: X %.3f Nmax %.2f vs N* %.2f -> %s; r %.3e' % (H, K, ns, c['X'], c['Nmax'], Nst,
                                                                                 'OK' if c['capacity_ok'] else 'FAILS', r))
    M.dump(M.RES / 'cmb_rescale.json', out)


if __name__ == '__main__':
    main(sys.argv[1:])

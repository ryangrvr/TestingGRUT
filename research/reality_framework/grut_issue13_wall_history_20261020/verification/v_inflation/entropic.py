"""Entropic (transverse-partition) factors R_X, R_D of the multi-field Langer rate (verifier v_inflation).

R = Z_perp(saddle phi = 0) / Z_perp(well phi = phi_m), Z_perp = int d^n z_perp exp(-(8 pi^2/(3 H^4)) V_perp),
with X complex = 2 real components (rho^2 = q1^2 + q2^2 = 2|X|^2, measure 2 pi rho drho) and D one real component.
V_X(phi, rho) = (1/2)(m2_X(mu) + g^2 phi^2 + 12 xi_X H^2) rho^2 + (lam_X/4) rho^4,
V_D(phi, D)   = (1/2)(m2_D(mu) + h^2 phi^2 + 12 xi_D H^2) D^2 + (lam_D/4) D^4, couplings at mu = H.
Also reports where the X fluctuation range at the saddle reaches the region in which phi = 0 stops being a saddle
of V(phi, |X|) (phi-curvature at phi=0: m2_phi + 12 xi_phi H^2 + g^2 rho^2 > 0 for rho^2 > rho_c^2).
"""
import json, math
import numpy as np
from scipy.integrate import quad
import rg_xi

V = 1000.0


def factors(H):
    r = rg_xi.run_to(np.array([max(H, rg_xi.mu0)]))
    c = rg_xi.couplings(*r[float(max(H, rg_xi.mu0))])
    S = 8 * math.pi ** 2 / (3 * H ** 4)
    a2 = c['m2_phi'] + 12 * c['xi_phi'] * H ** 2
    phim2 = -a2 / c['lam_phi']

    def ZX(phi2):
        m2 = c['m2_X'] + c['gX2'] * phi2 + 12 * c['xi_X'] * H ** 2
        f = lambda u: 0.5 * math.exp(-S * (0.5 * m2 * u + 0.25 * c['lam_X'] * u * u))  # u = rho^2, rho drho = du/2
        # scale of u
        us = min(1.0 / (S * abs(m2) + 1e-300), 1.0 / math.sqrt(S * c['lam_X']))
        return quad(f, 0, 200 * us, limit=400)[0]

    def ZD(phi2):
        m2 = c['m2_D'] + c['hD2'] * phi2 + 12 * c['xi_D'] * H ** 2
        f = lambda D: 2 * math.exp(-S * (0.5 * m2 * D * D + 0.25 * c['lam_D'] * D ** 4))
        ds = min(1.0 / math.sqrt(S * abs(m2)), (S * c['lam_D']) ** -0.25)
        return quad(f, 0, 60 * ds, limit=400)[0]
    RX = ZX(0.0) / ZX(phim2)
    RD = ZD(0.0) / ZD(phim2)
    # Gaussian well-only estimates
    RD_gauss = math.sqrt((c['m2_D'] + c['hD2'] * phim2) / c['m2_D'])
    # X fluctuation scale at the saddle vs critical rho^2
    rho_c2 = -a2 / c['gX2']
    # <rho^2> at saddle
    m2s = c['m2_X'] + 12 * c['xi_X'] * H ** 2
    num = quad(lambda u: u * math.exp(-S * (0.5 * m2s * u + 0.25 * c['lam_X'] * u * u)), 0, 200 / math.sqrt(S * c['lam_X']), limit=400)[0]
    den = quad(lambda u: math.exp(-S * (0.5 * m2s * u + 0.25 * c['lam_X'] * u * u)), 0, 200 / math.sqrt(S * c['lam_X']), limit=400)[0]
    rho2_mean = num / den
    pfrac = quad(lambda u: math.exp(-S * (0.5 * m2s * u + 0.25 * c['lam_X'] * u * u)), rho_c2, 200 / math.sqrt(S * c['lam_X']), limit=400)[0] / den
    # X mass at the well vs H
    mX_well = math.sqrt(c['m2_X'] + c['gX2'] * phim2 + 12 * c['xi_X'] * H ** 2)
    # X relaxation time at saddle (e-folds) ~ 3 H^2 / m_eff^2 with m_eff^2 = 3 lam_X <rho^2>
    tX = 3 * H ** 2 / (3 * c['lam_X'] * rho2_mean)
    return dict(H=H, R_X=RX, R_D=RD, R_D_gauss=RD_gauss, rho_c2_over_H2=rho_c2 / H ** 2,
                saddle_mean_rho2_over_H2=rho2_mean / H ** 2, P_saddle_rho2_gt_rho_c2=pfrac,
                mX_well_over_H=mX_well / H, X_relax_efolds_at_saddle=tX)


if __name__ == '__main__':
    out = [factors(H) for H in (100.0, 141.42, 251.19, 300.0, 335.0, 350.0, 500.0)]
    for o in out:
        print(o)
    json.dump(out, open('results_entropic.json', 'w'), indent=1)

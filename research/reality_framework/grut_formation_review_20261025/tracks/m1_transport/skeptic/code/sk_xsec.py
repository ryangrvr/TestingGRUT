"""Skeptic check of m1_transport T1-T3: tree-level cross sections recomputed by a DIFFERENT method.

Own code, independent of both Codex's elastic.py and the reviewer's xsec.py.  Instead of the
charge-basis Feynman rules (contact -4 lambda etc.) the amplitudes are built in the real-component
basis X = (phi1 + i phi2)/sqrt2 with ordinary real-scalar vertices (V_abcd = d^4 V/dphi^4,
h phi_a phi_b vertex 2 g^2 v delta_ab) and then rotated to the charge basis with e = (1, i)/sqrt2.
Units m_X = 1.  Potential (M1, d4747cb CANDIDATE_M1.json):
    V = lambda_phi/4 (phi^2 - v^2)^2 - eps v (phi - v)^3 + g^2 phi^2 |X|^2 + lambda_X |X|^4.
"""
import json, math
from pathlib import Path
import numpy as np

RES = Path(__file__).resolve().parents[1] / 'results'
g, lam, lam_phi, eps = 0.1, 1e-3, 0.03, 0.002
v_GeV = 0.924089637375362
mX = g * v_GeV
mphi2 = 2 * lam_phi * v_GeV**2 / mX**2            # 6 exactly; tilt has zero second derivative at v
g3 = 2 * g                                        # h phi_a phi_a vertex 2 g^2 v = 2 g m_X  (m_X units)

# real-basis quartic: lambda |X|^4 = lambda/4 (phi1^2 + phi2^2)^2
def V4(a, b, c, d):
    """fourth derivative of lambda/4 (phi1^2+phi2^2)^2 w.r.t. phi_a phi_b phi_c phi_d (indices 0/1)."""
    idx = sorted((a, b, c, d))
    n1 = idx.count(0)
    if n1 == 4 or n1 == 0:
        return 6 * lam          # lambda/4 * 4! = 6 lambda
    if n1 == 2:
        return 2 * lam          # lambda/2 * (2!)(2!) = 2 lambda
    return 0.0


def M_real(a, b, c, d, s, t, u):
    """tree amplitude phi_a phi_b -> phi_c phi_d (all real), m_X = 1."""
    m = -V4(a, b, c, d)
    m += g3**2 * ((a == b) * (c == d) / (mphi2 - s) + (a == c) * (b == d) / (mphi2 - t) + (a == d) * (b == c) / (mphi2 - u))
    return m


e = np.array([1, 1j]) / math.sqrt(2)      # particle X
eb = np.conj(e)                             # antiparticle Xbar


def amp(in1, in2, out1, out2, s, t, u):
    tot = 0j
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    tot += in1[a] * in2[b] * np.conj(out1[c]) * np.conj(out2[d]) * M_real(a, b, c, d, s, t, u)
    return tot


def main():
    out = {}
    s = 4.0
    # XX -> XX at threshold (t = u = 0)
    M_XX = amp(e, e, e, e, s, 0.0, 0.0)
    M_XX_contact_only = sum(e[a] * e[b] * np.conj(e[c]) * np.conj(e[d]) * (-V4(a, b, c, d))
                            for a in range(2) for b in range(2) for c in range(2) for d in range(2))
    sig_XX = abs(M_XX)**2 / (16 * math.pi * s) * 0.5          # identical final particles
    out['XX_to_XX_threshold'] = {
        'method': 'real-component basis (phi1, phi2) rotated to the charge basis',
        'M_total': float(M_XX.real), 'M_imag_residual': float(abs(M_XX.imag)),
        'M_contact_only': float(M_XX_contact_only.real),
        'M_h_exchange_t_plus_u': float((M_XX - M_XX_contact_only).real),
        'sigma_mX2': sig_XX, 'codex_sigma_mX2': 2.16627561e-7, 'reviewer_sigma_mX2': 2.16627561e-7,
        'ratio_to_codex': sig_XX / 2.16627561e-7,
        'sigma_contact_only_mX2': abs(M_XX_contact_only)**2 / (32 * math.pi * s),
        'note': 'Codex DERIVATIONS.md states the amplitude as 4 lambda_X - g_h^2/3, i.e. contact PLUS exchange; it never called sigma a contact cross section. The reviewer correction 1 addresses the task text, not Codex.'}
    # angular dependence at finite relative velocity (NR): t = -4p^2 (1-cos)/2 ... check threshold accuracy
    ang = {}
    for vrel in (0.0684, 0.2):
        p = vrel / 2
        s_v = 4 * (1 + p * p)
        cs = np.linspace(-1, 1, 2001)
        t = -2 * p * p * (1 - cs); u = -2 * p * p * (1 + cs)
        M = np.array([amp(e, e, e, e, s_v, tt, uu).real for tt, uu in zip(t, u)])
        sig = 0.5 * np.trapezoid(M**2, cs) / (32 * math.pi * s_v)   # dsigma/dcos = |M|^2/(32 pi s) (before 1/2)
        ang[str(vrel)] = {'sigma_mX2': float(sig), 'ratio_to_threshold': float(sig / sig_XX)}
    out['XX_velocity_dependence'] = ang
    # X Xbar -> X Xbar at threshold: in X(a) Xbar(b) -> X(c) Xbar(d); t = (p_a - p_c)^2 -> 0, u -> 0
    M_XXb = amp(e, eb, e, eb, s, 0.0, 0.0)
    out['XXbar_to_XXbar_threshold'] = {'M_total': float(M_XXb.real), 'M_imag_residual': float(abs(M_XXb.imag)),
                                       'sigma_mX2': abs(M_XXb)**2 / (16 * math.pi * s),
                                       'reviewer_sigma_mX2': 2.555e-6,
                                       's_channel_term': g3**2 / (mphi2 - s), 't_channel_term': g3**2 / mphi2,
                                       'contact_term': float(sum(e[a] * eb[b] * np.conj(e[c]) * np.conj(eb[d]) * (-V4(a, b, c, d))
                                                                 for a in range(2) for b in range(2) for c in range(2) for d in range(2)).real)}
    # Gamma(h -> X Xbar) from the real basis: h -> phi_a phi_a, identical real daughters, summed over a
    beta = math.sqrt(1 - 4 / mphi2)
    Gam_real = 2 * g3**2 * beta / (32 * math.pi * math.sqrt(mphi2))     # 2 channels x identical-particle width
    out['Gamma_h_to_XXbar'] = {'mX_units': Gam_real, 'GeV': Gam_real * mX, 'codex_mX_units': 1.875659e-4,
                               'reviewer_mX_units': 1.875659e-4, 'ratio': Gam_real / 1.875659e-4,
                               'm_phi_over_mX': math.sqrt(mphi2), 'tilt_second_derivative_at_v': 0.0}
    # NR Yukawa from h exchange: V(r) = -(g3^2/4m^2) e^{-m_h r}/(4 pi r) = -alpha e^{-m_h r}/r
    alpha = g3**2 / (4 * 1.0) / (4 * math.pi)
    out['NR_potential_checks'] = {'alpha_yukawa': alpha, 'alpha_mX_over_mphi': alpha / math.sqrt(mphi2),
                                  'bound_state_threshold_alpha_m_over_mphi': 0.84,
                                  'sommerfeld_alpha_over_v_at_0.0684': alpha / 0.0684,
                                  'born_valid': True}
    out['SIDM'] = {'sigma_cm2': sig_XX / mX**2 * 0.3894e-27, 'sigma_over_m_cm2_g': sig_XX / mX**2 * 0.3894e-27 / (mX * 1.78266e-24)}
    RES.mkdir(parents=True, exist_ok=True)
    (RES / 'sk_xsec.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()

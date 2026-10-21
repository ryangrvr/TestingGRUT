"""Independent one-loop RG of the Candidate 1 scalar sector including nonminimal couplings.

Verifier: v_inflation.  Written from scratch (no module code imported or copied).

Method: four real fields z = (phi, x, y, D), X = (x + i y)/sqrt2.  The full symmetric quartic tensor
Lam_abcd = d^4 V / dz_a dz_b dz_c dz_d, the mass matrix M_ab (V2 = M_ab z_a z_b / 2) and the nonminimal
matrix Xi_ab (V_R = Xi_ab R z_a z_b / 2, conformal value Xi = delta/6) are integrated directly as tensors:

  16 pi^2 dLam_abcd/dt = Lam_abef Lam_efcd + Lam_acef Lam_efbd + Lam_adef Lam_efbc
  16 pi^2 dM_ab/dt    = Lam_abcd M_cd
  16 pi^2 dXi_ab/dt   = Lam_abcd (Xi_cd - delta_cd/6)

t = ln mu.  One-loop scalar anomalous dimension is zero (pure quartic scalar theory, no gauge/Yukawa).
The xi formula is the standard curved-space result (Buchbinder-Odintsov-Shapiro book, not read here;
normalization checked against Herranen-Markkanen-Nurmi-Rajantie 1407.3141 eq. for the SM Higgs,
16 pi^2 beta_xi = (xi - 1/6)(12 lambda + ...), for V = lambda (H^dag H)^2: reproduced below as a test).

The quartic tensor is built from the polynomial by exact symbolic differentiation of monomials
(no finite differences).
"""
import itertools, json, math
import numpy as np
from scipy.integrate import solve_ivp

NF = 4
PI2_16 = 16.0 * math.pi ** 2


def quartic_tensor(monos):
    """monos: dict {exponent tuple (len 4): coefficient}. returns Lam_abcd = 4th derivative tensor."""
    L = np.zeros((NF,) * 4)
    for idx in itertools.product(range(NF), repeat=4):
        cnt = [0] * NF
        for i in idx:
            cnt[i] += 1
        key = tuple(cnt)
        if key in monos:
            # d^4/dz^cnt of coef * prod z^cnt = coef * prod cnt_i!
            L[idx] = monos[key] * np.prod([math.factorial(c) for c in cnt])
    return L


def candidate_monomials(lp, lx, ld, c, d, alpha):
    """V4 = lp/4 phi^4 + lx |X|^4 + ld/4 D^4 + c phi^2|X|^2 + d/2 phi^2 D^2 + alpha |X|^2 D^2,
    with |X|^2 = (x^2+y^2)/2.  (c = g_X^2, d = h_D^2.)"""
    m = {}
    def add(k, v):
        m[k] = m.get(k, 0.0) + v
    add((4, 0, 0, 0), lp / 4)
    # lx |X|^4 = lx (x^2+y^2)^2 / 4
    add((0, 4, 0, 0), lx / 4); add((0, 0, 4, 0), lx / 4); add((0, 2, 2, 0), lx / 2)
    add((0, 0, 0, 4), ld / 4)
    # c phi^2 |X|^2 = c/2 phi^2 (x^2 + y^2)
    add((2, 2, 0, 0), c / 2); add((2, 0, 2, 0), c / 2)
    add((2, 0, 0, 2), d / 2)
    # alpha |X|^2 D^2 = alpha/2 (x^2+y^2) D^2
    add((0, 2, 0, 2), alpha / 2); add((0, 0, 2, 2), alpha / 2)
    return m


def beta_lam(L):
    t1 = np.einsum('abef,efcd->abcd', L, L)
    t2 = np.einsum('acef,efbd->abcd', L, L)
    t3 = np.einsum('adef,efbc->abcd', L, L)
    return (t1 + t2 + t3) / PI2_16


def beta_M(L, M):
    return np.einsum('abcd,cd->ab', L, M) / PI2_16


def beta_xi(L, Xi):
    return np.einsum('abcd,cd->ab', L, Xi - np.eye(NF) / 6.0) / PI2_16


def pack(L, M, Xi):
    return np.concatenate([L.ravel(), M.ravel(), Xi.ravel()])


def unpack(y):
    L = y[:256].reshape((NF,) * 4)
    M = y[256:272].reshape(NF, NF)
    Xi = y[272:288].reshape(NF, NF)
    return L, M, Xi


def rhs(t, y):
    L, M, Xi = unpack(y)
    return pack(beta_lam(L), beta_M(L, M), beta_xi(L, Xi))


# benchmark (BRIEF.md)
v = 1000.0
gX, lp, lx, hD, ld, mD0 = 0.1, 0.01, 0.001, 0.01, 0.001, 100.0
mu0 = gX * v  # 100 GeV


def initial():
    L0 = quartic_tensor(candidate_monomials(lp, lx, ld, gX ** 2, hD ** 2, 0.0))
    M0 = np.diag([-lp * v ** 2, 0.0, 0.0, mD0 ** 2])  # m_X0^2 = 0 at mu0
    Xi0 = np.zeros((NF, NF))
    return L0, M0, Xi0


def run_to(mus, rtol=1e-12, atol=1e-16):
    L0, M0, Xi0 = initial()
    ts = np.log(np.asarray(mus) / mu0)
    if ts.max() <= 0.0:
        return {float(m): (L0.copy(), M0.copy(), Xi0.copy()) for m in mus}
    order = np.argsort(ts)
    sol = solve_ivp(rhs, (0.0, max(ts.max(), 0.0)), pack(L0, M0, Xi0), t_eval=np.sort(ts),
                    method='DOP853', rtol=rtol, atol=atol)
    out = {}
    for j, k in enumerate(order):
        out[float(mus[k])] = unpack(sol.y[:, j])
    return out


def couplings(L, M, Xi):
    return dict(lam_phi=L[0, 0, 0, 0] / 6, lam_X=L[1, 1, 1, 1] / 6, lam_D=L[3, 3, 3, 3] / 6,
                gX2=L[0, 0, 1, 1] / 2, hD2=L[0, 0, 3, 3] / 2, alpha=L[1, 1, 3, 3] / 2,
                m2_phi=M[0, 0], m2_X=M[1, 1], m2_D=M[3, 3],
                xi_phi=Xi[0, 0], xi_X=Xi[1, 1], xi_D=Xi[3, 3])


def tests():
    res = {}
    # SM-Higgs-like O(4) test of the xi normalization: V = lam (h.h/2)^2 = lam/4 (sum z^2)^2
    lam = 0.13
    m = {}
    for i in range(4):
        e = [0] * 4; e[i] = 4; m[tuple(e)] = lam / 4
    for i, j in itertools.combinations(range(4), 2):
        e = [0] * 4; e[i] = 2; e[j] = 2; m[tuple(e)] = lam / 2
    L = quartic_tensor(m)
    xi = 0.3
    b = beta_xi(L, xi * np.eye(4))[0, 0] * PI2_16
    res['O4_beta_xi_over_(xi-1/6)'] = b / (xi - 1 / 6)
    res['O4_expected_12lambda'] = 12 * lam
    # O(4) quartic beta: known 16pi^2 beta_lambda = 24 lambda^2 for V = lambda (H^dag H)^2 (pure scalar)
    bl = beta_lam(L)[0, 0, 0, 0] * PI2_16 / 6
    res['O4_beta_lambda_over_lambda2'] = bl / lam ** 2
    res['O4_expected_24'] = 24.0
    # Candidate closed forms at mu0 (derived by hand in VERIFY.md)
    L0, M0, Xi0 = initial()
    B = beta_lam(L0) * PI2_16
    res['B_lp_tensor'] = B[0, 0, 0, 0] / 6
    res['B_lp_hand'] = 18 * lp ** 2 + 4 * gX ** 4 + 2 * hD ** 4
    res['B_alpha_tensor'] = B[1, 1, 3, 3] / 2
    res['beta_alpha_mu0'] = B[1, 1, 3, 3] / 2 / PI2_16
    bx = beta_xi(L0, Xi0) * PI2_16
    res['B_xi_phi_mu0_tensor'] = bx[0, 0]
    res['B_xi_phi_mu0_hand'] = -(6 * lp + 4 * gX ** 2 + 2 * hD ** 2) / 6
    res['B_xi_X_mu0_tensor'] = bx[1, 1]
    res['B_xi_D_mu0_tensor'] = bx[3, 3]
    res['beta_xi_phi_mu0'] = bx[0, 0] / PI2_16
    res['B_m2phi_mu0'] = (beta_M(L0, M0) * PI2_16)[0, 0]
    return res


if __name__ == '__main__':
    out = {'tests': tests()}
    Hs = [1e2, 2.5e2, 3.5e2, 1e3, 2e3, 1e4, 1e6, 1e8, 1e10, 1e13]
    runs = run_to(np.array(Hs))
    tab = {}
    for H in Hs:
        cpl = couplings(*runs[float(H)])
        cpl['12xi_phi'] = 12 * cpl['xi_phi']
        tab['%.4g' % H] = cpl
    out['running_mu_eq_H'] = tab
    # leading log
    out['leading_log_xi_phi_1e13'] = out['tests']['beta_xi_phi_mu0'] * math.log(1e13 / mu0)
    # scale band at 1e13
    band = run_to(np.array([0.5e13, math.sqrt(12) * 1e13]))
    out['xi_phi_1e13_band_H/2_sqrt12H'] = [couplings(*band[k])['xi_phi'] for k in sorted(band)]
    # convention conversion: action text writes xi R phi^2 (no 1/2)
    out['xi_phi_1e13_in_action_text_convention_(xi R phi^2)'] = tab['1e+13']['xi_phi'] / 2
    # phi_m / H at 1e13 on X=0 slice
    c13 = tab['1e+13']
    a = -(c13['m2_phi'] / 1e26 + 12 * c13['xi_phi'])
    out['phi_m_over_H_1e13'] = math.sqrt(a / c13['lam_phi'])
    json.dump(out, open('results_rg.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))

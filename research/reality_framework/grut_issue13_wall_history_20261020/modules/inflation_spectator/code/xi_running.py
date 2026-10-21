"""Part 1: one-loop running of the nonminimal couplings from xi(mu_0) = 0 (matching) to mu ~ H_inf.

Derivation (NOTES.md section 1): for real fields z_a with V4 = Lambda_abcd z_a z_b z_c z_d / 4! and
nonminimal term (1/2) xi_ab R z_a z_b, the one-loop curvature-linear divergence is obtained from the
mass-squared one by m^2_ab -> (xi_ab - delta_ab/6) R, giving
    16 pi^2 beta_{xi_ab} = sum_cd Lambda_abcd (xi_cd - delta_cd/6).
Checks performed here:
  (a) closed forms in common.B_all vs the tensor contraction, with Lambda from exact finite-difference
      fourth derivatives of the quartic polynomial (independent of the hand-written entries);
  (b) the curvature-linear RG identity sum_ab beta_xi_ab z_a z_b /2 = Tr[H4(z) Xi]/(16 pi^2), H4 = quartic part of the Hessian, at random
      fields and random couplings;
  (c) single O(4) field check against Herranen-Markkanen-Nurmi-Rajantie arXiv:1407.3141 eq. (7), whose
      scalar part is 16 pi^2 beta_xi = 12 lambda (xi - 1/6) for V = lambda (H^dag H)^2.
"""
import itertools
import numpy as np
from common import B_all, couplings, K16, MU0, NAMES, P0, dump


def quartic(z, lp, lx, ld, c, d, a):
    ph, x, y, D = z
    R2 = x * x + y * y
    return (lp * ph ** 4 / 4 + lx * R2 ** 2 / 4 + ld * D ** 4 / 4 + c * ph ** 2 * R2 / 2
            + d * ph ** 2 * D ** 2 / 2 + a * R2 * D ** 2 / 2)


def lambda_tensor(f, n=4, h=0.5):
    """Fourth derivatives of a quartic polynomial by the exact mixed central difference
    d^4 f / dz_a dz_b dz_c dz_d = sum over 16 sign choices of s_a s_b s_c s_d f(sum s_i h e_i) / (2h)^4,
    which is exact (up to rounding) for a homogeneous quartic."""
    L = np.zeros((n, n, n, n))
    E = np.eye(n)
    for idx in itertools.product(range(n), repeat=4):
        if L[idx] != 0:
            continue
        tot = 0.0
        for sg in itertools.product((1, -1), repeat=4):
            pt = sum(sg[k] * h * E[idx[k]] for k in range(4))
            tot += np.prod(sg) * f(pt)
        val = tot / (2 * h) ** 4
        for perm in set(itertools.permutations(idx)):
            L[perm] = val
    return L


def checks():
    rng = np.random.default_rng(20261021)
    out = {}
    # (a) tensor vs closed form, at the matching point and at random couplings / xis
    worst = 0.0
    for trial in range(6):
        if trial == 0:
            p = P0.copy()
        else:
            p = P0.copy()
            p[:6] = rng.uniform(0.0005, 0.05, 6)
            p[9:] = rng.uniform(-0.3, 0.3, 3)
        lp, lx, ld, c, d, a = p[:6]
        L = lambda_tensor(lambda z: quartic(z, lp, lx, ld, c, d, a))
        xi = np.diag([p[9], p[10], p[10], p[11]]) - np.eye(4) / 6.0
        Bt = np.einsum("abcd,cd->ab", L, xi)
        B = B_all(p)
        closed = np.diag([B[9], B[10], B[10], B[11]])
        offdiag = np.max(np.abs(Bt - np.diag(np.diag(Bt))))
        dev = np.max(np.abs(Bt - closed)) / max(1e-30, np.max(np.abs(closed)))
        worst = max(worst, dev, offdiag)
    out["tensor_vs_closed_form_max_rel_dev"] = worst
    # (b) curvature-linear RG identity at random fields
    worst_b = 0.0
    for trial in range(20):
        p = P0.copy()
        p[:6] = rng.uniform(0.0005, 0.05, 6)
        p[9:] = rng.uniform(-0.3, 0.3, 3)
        lp, lx, ld, c, d, a = p[:6]
        z = rng.normal(size=4)
        ph, x, y, D = z
        R2 = x * x + y * y
        # quartic-generated part of the Hessian H4_ab = d^2 V4 / dz_a dz_b
        H4 = np.array([
            [3 * lp * ph ** 2 + c * R2 + d * D ** 2, 2 * c * ph * x, 2 * c * ph * y, 2 * d * ph * D],
            [2 * c * ph * x, c * ph ** 2 + lx * (3 * x * x + y * y) + a * D ** 2, 2 * lx * x * y, 2 * a * x * D],
            [2 * c * ph * y, 2 * lx * x * y, c * ph ** 2 + lx * (x * x + 3 * y * y) + a * D ** 2, 2 * a * y * D],
            [2 * d * ph * D, 2 * a * x * D, 2 * a * y * D, d * ph ** 2 + a * R2 + 3 * ld * D ** 2]])
        Xi = np.diag([p[9], p[10], p[10], p[11]]) - np.eye(4) / 6.0
        rhs = 2.0 * np.trace(H4 @ Xi)          # sum_ab B_ab z_a z_b = sum Lambda_abcd z_a z_b Xi_cd = 2 Tr[H4 Xi]
        B = B_all(p)
        lhs = B[9] * ph ** 2 + B[10] * R2 + B[11] * D ** 2      # = sum_ab B_xi_ab z_a z_b
        worst_b = max(worst_b, abs(lhs - rhs) / max(abs(rhs), 1e-30))
    out["curvature_linear_RG_identity_max_rel_dev"] = worst_b
    # (c) single O(4) multiplet with V = lam (h^2/2)^2: B_xi = (xi-1/6) * (Lambda_1111 + 3 Lambda_1122)
    lam = 0.13
    L = lambda_tensor(lambda z: lam * (0.5 * np.dot(z, z)) ** 2)
    coeff = L[0, 0, 0, 0] + L[0, 0, 1, 1] + L[0, 0, 2, 2] + L[0, 0, 3, 3]
    out["O4_check_coefficient_over_lambda"] = coeff / lam
    out["O4_check_expected_HMNR_eq7"] = 12.0
    return out


def main():
    out = {"module": "inflation_spectator", "part": "1 xi running",
           "convention": "potential term +(1/2) xi_phi R phi^2 + xi_X R |X|^2 + (1/2) xi_D R D^2; R = 12 H^2 (de Sitter, -+++); "
                         "Hubble mass^2 = 12 xi H^2; xi(mu_0 = 100 GeV) = 0 by matching (GRUT_CANDIDATE_1.md section 3)",
           "label": "computed (one-loop pure-scalar RG; SM does not enter because all portals vanish at mu_0 and are "
                    "multiplicatively renormalised at one loop; graviton-loop terms ~ (H/M_P)^2 not included)"}
    out["checks"] = checks()
    B0 = B_all(P0)
    out["B_xi_at_mu0"] = {"xi_phi": B0[9], "xi_X": B0[10], "xi_D": B0[11]}
    out["beta_xi_at_mu0"] = {"xi_phi": B0[9] / K16, "xi_X": B0[10] / K16, "xi_D": B0[11] / K16}
    out["leading_log_formula"] = ("xi_phi(mu) ~= -(6 lp + 4 c + 2 d)/(6 * 16 pi^2) * ln(mu/mu_0) "
                                  "= %.6e * ln(mu/mu_0)" % (B0[9] / K16))
    rows = []
    for lg in np.arange(2.0, 13.01, 0.5):
        H = 10 ** lg
        row = {"H_GeV": H}
        for tag, k in (("muH", 1.0), ("muSqrtR", np.sqrt(12.0)), ("muHhalf", 0.5)):
            if k * H < MU0 * np.exp(-3.0):
                continue
            p = couplings(k * H)
            row[tag] = {"xi_phi": p["xp"], "xi_X": p["xx"], "xi_D": p["xd"],
                        "hubble_mass2_over_H2_phi": 12 * p["xp"], "hubble_mass2_over_H2_X": 12 * p["xx"],
                        "hubble_mass2_over_H2_D": 12 * p["xd"],
                        "lam_phi": p["lp"], "lam_X": p["lx"], "gX2": p["c"], "mphi2_param_GeV2": p["mp"],
                        "mX0_2_generated_GeV2": p["mx"], "mD2_param_GeV2": p["md"]}
        row["leading_log_xi_phi"] = B0[9] / K16 * np.log(H / MU0)
        rows.append(row)
    out["table"] = rows
    pH = couplings(1e13)
    out["headline"] = {
        "xi_phi_range_muH_H_1e2_to_1e13": [couplings(1e2)["xp"], pH["xp"]],
        "xi_phi_at_H_1e13_scale_band_[H/2,H,sqrt(12)H]": [couplings(0.5e13)["xp"], pH["xp"], couplings(np.sqrt(12) * 1e13)["xp"]],
        "hubble_mass2_over_H2_phi_at_1e13": 12 * pH["xp"],
        "sign": "xi_phi(mu) < 0 for mu > mu_0: the curvature-induced phi mass^2 is tachyonic (adds to the double well); "
                "|12 xi_phi| <= 0.033, so phi stays light (|m_H^2| << H^2) and the shift is small",
        "D_tachyonic_above_H_GeV": None,
    }
    # D: m_D^2(phi=0) = m_D0^2 + 12 xi_D H^2 turns negative above some H
    Hs = np.logspace(2, 13, 2000)
    md_eff = np.array([couplings(h)["md"] + 12 * couplings(h)["xd"] * h ** 2 for h in Hs])
    neg = Hs[md_eff < 0]
    out["headline"]["D_tachyonic_above_H_GeV"] = float(neg[0]) if neg.size else None
    mx_eff = np.array([couplings(h)["mx"] + 12 * couplings(h)["xx"] * h ** 2 for h in Hs])
    negx = Hs[mx_eff < 0]
    out["headline"]["X_mass2_at_phi0_negative_above_H_GeV"] = float(negx[0]) if negx.size else None
    dump("xi_running.json", out)
    return out


if __name__ == "__main__":
    o = main()
    print(o["checks"])
    print(o["B_xi_at_mu0"], o["beta_xi_at_mu0"])
    print(o["headline"])
    for r in o["table"][::4]:
        print(r["H_GeV"], r["muH"]["xi_phi"], r["muH"]["xi_X"], r["muH"]["xi_D"], r["muH"]["lam_phi"])

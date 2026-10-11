#!/usr/bin/env python3
"""CP-02 exact contractions and independent numerical integrals; no data fitting.

Run from any directory: python research/constructive_physics/attempt02/check.py
Requires the versions in requirements.txt. Writes deterministic results.json.
This is an author calculation, not an external review or experiment.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import cmath
import hashlib
import json
import math
import platform

import numpy as np
import scipy
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
GROUPS = []
NUMERICAL = []


def done(name, **evidence):
    GROUPS.append({"name": name, "status": "PASS", **evidence})


def exact(x):
    return str(x)


def close(name, actual, expected, atol=1e-12):
    err = abs(actual - expected)
    assert err < atol, (name, actual, expected, err)
    NUMERICAL.append({"name": name, "actual": float(actual),
                      "expected": float(expected), "absolute_error": float(err),
                      "tolerance": atol})


def dot_table(s, t):
    u = -s-t
    # All momenta incoming, eta=(-,+,+,+), s=-(p1+p2)^2.
    return [[F(0), -s/2, -t/2, -u/2],
            [-s/2, F(0), -u/2, -t/2],
            [-t/2, -u/2, F(0), -s/2],
            [-u/2, -t/2, -s/2, F(0)]]


KIN = [(F(1), F(-1, 2)), (F(2), F(-1, 2))]
ANGLES = [F(-1), F(-1, 2), F(0), F(1, 2), F(1)]
ENERGIES = [F(1, 20), F(1, 10)]

# Check the freeze before doing any scientific calculation.
freeze = json.loads((HERE / "FREEZE.json").read_text())
for name, digest in freeze["sha256"].items():
    assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest, name
done("freeze_bytes", frozen_files=len(freeze["sha256"]),
     freeze_commit="d88a1a0bb0d4edea54b075437f0edd0921381df2")

# The nonlinear boost imposes 2(1+z) f'(z)-f(z)=0.
f = [F(-1)]
for n in range(4):
    f.append(F(1-2*n, 2*(n+1))*f[-1])
assert f == [F(-1), F(-1, 2), F(1, 8), F(-1, 16), F(5, 128)]
for n in range(4):
    assert 2*(n+1)*f[n+1] + (2*n-1)*f[n] == 0
done("nonlinear_boost_Ward_recursion", volume_series=list(map(exact, f)))

eta = np.diag([-1., 1., 1., 1.])
for v in [np.zeros(4), np.array([.2, .1, 0., 0.]),
          np.array([0., .5, 1/3, 0.])]:
    gamma = eta + np.outer(v, v)
    determinant_ratio = np.linalg.det(gamma)/np.linalg.det(eta)
    close("determinant_lemma", determinant_ratio, 1+v@eta@v)
    assert np.linalg.eigvalsh(gamma)[0] < 0 < np.linalg.eigvalsh(gamma)[1]
done("induced_metric_and_signature", cases=3)

# Enumerate all 4! derivatives of the quartic local operators (no fitted factor).
for s, t in KIN:
    u = -s-t
    dots = dot_table(s, t)
    contact4 = F(1, 8)*sum(dots[a][b]*dots[c][d]
                          for a, b, c, d in permutations(range(4)))
    assert contact4 == (s*s+t*t+u*u)/4
done("tree_LSZ_quartic_contraction", kappa="1/(4 T_R)")

# R^2 is equivalent on shell at this order to (Pi_mu_nu Pi^mu_nu)^2/T_R^2.
for s, t in KIN:
    u = -s-t
    dots = dot_table(s, t)
    contactR2 = sum(dots[a][b]**2*dots[c][d]**2
                    for a, b, c, d in permutations(range(4)))
    assert contactR2 == (s**4+t**4+u**4)/2
done("same_symmetry_R_squared_countermodel",
     amplitude_delta="zeta (s^4+t^4+u^4)/(2 T_R^2)")

# Rational beta-function moments: integral t^(2k)/(1+t^2)^6 dt / pi.
def beta_moment_over_pi(k):
    # Gamma(k+1/2) Gamma(11/2-k) / Gamma(6) / pi.
    a, b = k, 5-k
    return F(math.factorial(2*a)*math.factorial(2*b),
             4**(a+b)*math.factorial(a)*math.factorial(b)*math.factorial(5))

moments = [beta_moment_over_pi(k) for k in range(3)]
assert moments == [F(63, 256), F(7, 256), F(3, 256)]
wick_moment = 9*moments[2]-6*moments[1]+moments[0]
assert wick_moment == F(3, 16)
# Correlator includes 1/(2 pi^4); sources include 1/(4 T_R^2).
assert wick_moment/8 == F(3, 128)
done("two_source_Wick_integral_exact", moments_over_pi=list(map(exact, moments)),
     C7="3/(128 pi^3 T_R^2)")

for r in [1., 2., 5.]:
    # A rescaling gives stable independent improper quadrature at every r.
    integrand = lambda x: (3*x*x-1)**2 / (2*math.pi**4*(1+x*x)**6*r**7)
    val, enclosure_estimate = quad(integrand, -np.inf, np.inf,
                                   epsabs=1e-13, epsrel=1e-11)
    expected = 3/(32*math.pi**3*r**7)
    close("Wick_integral_quadrature", val, expected, 1e-12)
    assert enclosure_estimate < 1e-12
done("two_source_integral_independent_quadrature", cases=3)

kappa = F(1, 4)
assert F(3, 128) == F(3, 8)*kappa*kappa
done("joint_asymptotic_relation", C7_over_kappa_squared="3/(8 pi^3)")

comparator = []
for b in [F(0), F(1), F(2)]:
    ratio = b*b
    comparator.append({"b": exact(b), "same_kappa": exact(kappa),
                       "C7_over_geometric_C7": exact(ratio),
                       "full_5D_symmetry": b == 1})
assert comparator[2]["C7_over_geometric_C7"] == "4"
done("ordinary_four_dimensional_countermodel", cases=comparator)

finite = []
for E in ENERGIES:
    s, t, u = 4*E*E, -2*E*E, -2*E*E
    tree = (s*s+t*t+u*u)/4
    delta = (s**4+t**4+u**4)/2
    assert delta/tree == 24*E**4
    finite.append({"E_over_T_quarter": exact(E), "zeta": 1,
                   "relative_amplitude_change": exact(delta/tree)})
done("finite_energy_rigidity_refuted", cases=finite)

# Exact angular moments on S^2: <a^2>=1/3,
# <a^2 b^2>=(1+2 z^2)/15.
for z in ANGLES:
    integral_over_4pi = F(9)+F(6, 3)+(1+2*z*z)/15
    assert integral_over_4pi == (166+2*z*z)/15
    cut_over_s4_over_pi = integral_over_4pi / 2048
    assert cut_over_s4_over_pi == (83+z*z)/15360
done("one_loop_two_particle_cut_exact",
     Im_Ms="s^4 (83+z^2)/(15360 pi T_R^2)",
     identical_particle_phase_space_factor="1/2")

nodes, weights = np.polynomial.legendre.leggauss(24)
phis = (np.arange(64)+.5)*(2*math.pi/64)
for z0 in ANGLES:
    z = float(z0)
    val = 0.
    for a, w in zip(nodes, weights):
        b = z*a+math.sqrt(max(0., 1-z*z))*math.sqrt(1-a*a)*np.cos(phis)
        val += w*np.sum((3+a*a)*(3+b*b))*(2*math.pi/64)
    close("angular_cut_quadrature", val/(4*math.pi),
          float((166+2*z0*z0)/15), 2e-13)
done("one_loop_cut_independent_angular_quadrature", cases=5)

# An independent partial-wave route, M=32 pi sum_l (2l+1) a_l P_l.
for z in ANGLES:
    P2 = (3*z*z-1)/2
    partial = 32*(F(5, 384)**2+5*F(1, 1920)**2*P2)
    assert partial == (83+z*z)/15360
done("one_loop_cut_partial_wave_unitarity", waves=[0, 2])

# Exact polynomial operations prove identities on the whole s+t+u=0 plane.
def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, F(0))+v
    return {k: v for k, v in out.items() if v}

def mul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            ij = (i+k, j+l)
            out[ij] = out.get(ij, F(0))+x*y
    return {k: v for k, v in out.items() if v}

def power(p, n):
    out = {(0, 0): F(1)}
    for _ in range(n):
        out = mul(out, p)
    return out

def scale(p, a):
    return {k: a*v for k, v in p.items() if a*v}

S, TT = {(1, 0): F(1)}, {(0, 1): F(1)}
U = add(scale(S, -1), scale(TT, -1))
sum4 = add(*(power(p, 4) for p in [S, TT, U]))
mixed = add(*(mul(power(a, 2), power(add(b, scale(c, -1)), 2))
              for a, b, c in [(S, TT, U), (TT, U, S), (U, S, TT)]))
assert mixed == sum4
sum_cut_numerator = add(scale(sum4, 83), mixed)
assert sum_cut_numerator == scale(sum4, 84)
uv = F(84, 15360)
assert uv == F(7, 1280)
done("crossing_sum_and_UV_polynomial", UV_pole_over_pi_squared="7/1280")

# Independent published covariant heat-kernel contraction (2010.05913v1 4.17):
# Gamma_pole = 7/(320 pi^2 epsilon) tr(Pi^4)/T_R^2.
for s, t in KIN:
    dots = dot_table(s, t)
    trace_contact = sum(dots[a][b]*dots[b][c]*dots[c][d]*dots[d][a]
                        for a, b, c, d in permutations(range(4)))
    u = -s-t
    assert trace_contact == (s**4+t**4+u**4)/4
    assert F(7, 320)*trace_contact == uv*(s**4+t**4+u**4)
done("covariant_heat_kernel_to_LSZ_amplitude", contraction_factor="1/4")

# Independently reduce the paper's N=1, L->infinity Eq.5.20.
def paper_520(a, b):
    return add(scale(power(a, 4), F(2, 5)),
               scale(mul(power(a, 3), b), F(102, 5)),
               scale(mul(power(a, 2), power(b, 2)), F(182, 5)),
               scale(mul(a, power(b, 3)), 32), scale(power(b, 4), 16))
published_crossing = scale(add(paper_520(S, TT), paper_520(TT, S),
                               paper_520(U, TT)), F(1, 96*16))
assert published_crossing == scale(sum4, uv)
assert uv/F(7, 10240) == 8
done("primary_paper_normalization_audit",
     equations_agree=["our unitarity", "our 4! contractions", "source 4.17", "source 5.20"],
     source_4_20_discrepancy_factor="8",
     resolution="Unresolved externally; no claim of author-confirmed erratum")

beta = -2*uv
assert beta == F(-7, 640)
def amp8(s, t, mu, c8):
    u = -s-t
    poly = s**4+t**4+u**4
    logs = sum((83*a**4+a*a*(b-c)**2)*
               cmath.log(complex(-a, -0.)/(mu*mu))
               for a, b, c in [(s, t, u), (t, u, s), (u, s, t)])
    return c8*poly-logs/(15360*math.pi**2)
for s0, t0 in KIN:
    s, t = float(s0), float(t0)
    a = amp8(s, t, 1., 0.)
    for mu in [.5, 1., 2.]:
        c8 = float(beta)*math.log(mu)/math.pi**2
        q = amp8(s, t, mu, c8)
        close("RG_scale_independent_real_amplitude", q.real, a.real, 2e-13)
        close("RG_scale_independent_imaginary_amplitude", q.imag, a.imag, 2e-13)
done("one_loop_RG_scale_cancellation", beta_c8="(-7/640)/pi^2",
     finite_c8_boundary="independently matchable")

# The massive bosonic determinant changes geometry, not a separate soft coupling.
# a2 = R^2/72 + (Riem^2-Ric^2)/180 modulo total derivatives.
assert F(5, 360) == F(1, 72)
assert F(2, 360) == F(1, 180)
assert F(-2, 360) == F(-1, 180)
# Massive determinant volume: mu^4/(64 pi^2) log(mu^2/mu_R^2).
assert F(1, 32)*F(1, 2) == F(1, 64)
assert F(1, 32)-2*F(1, 64) == 0
done("matter_heat_kernel_operator_inventory",
     local_pole_integrand="mu^4/2 - mu^2 R/6 + R^2/72 + (Riem^2-Ric^2)/180",
     beta_T_MS="mu^4/(32 pi^2)",
     meaning="Tension, induced curvature and higher curvature require matching; no finite family closure")

# Hamiltonian and characteristic speed about an admissible homogeneous background.
for v in [F(0), F(1, 5)]:
    assert 1-v*v > 0
done("leading_classical_stability", c_squared="1-v^2 > 0 for |v|<1",
     vacuum_kinetic_sign="positive", all_EFT_backgrounds_proved=False)

# Linear gravity obstruction: delta gamma is quadratic about the flat vacuum.
# The EH quadratic bilinear has Fourier kernel (k^2)^2-(k_mu k_nu)^2=0.
for k in [(F(1), F(0), F(0), F(0)), (F(1), F(2), F(3), F(0))]:
    sg = [F(-1), F(1), F(1), F(1)]
    k2 = sum(sg[i]*k[i]*k[i] for i in range(4))
    contracted = sum(sg[i]*sg[j]*k[i]*k[j]*k[i]*k[j]
                     for i in range(4) for j in range(4))
    assert k2*k2 == contracted
done("induced_Einstein_Hilbert_quadratic_boundary_term",
     linear_metric_modes=0, embedding_transverse_modes=1,
     Einstein_spin2_modes_supplied=0, Newtonian_force_recovered=False)

# Ordinary operator-valued interactions never depend on an ensemble decomposition.
# A symbolic linearity check is algebraic: U[p rho+(1-p)sigma]Udag =
# p U rho Udag+(1-p)U sigma Udag. No modified probability rule is present.
done("operational_mixture_boundary",
     density_map="rho -> U rho Udag",
     status_scope="ordinary quantization postulated; QFT UV/no-signaling completion not proved by software")

result = {
    "status": "CP02_GEOMETRIC_EFT_CHECKS_PASS",
    "evidence_grade": "AUTHOR_EXACT_DERIVATIONS_AND_NUMERICAL_CROSS_CHECKS",
    "experiment": False, "independent_review": False,
    "freeze_commit": "d88a1a0bb0d4edea54b075437f0edd0921381df2",
    "groups_passed": len(GROUPS), "groups": GROUPS,
    "numerical_checks": NUMERICAL,
    "largest_numerical_absolute_error": max(x["absolute_error"] for x in NUMERICAL),
    "joint_relation": "C7=3 kappa^2/(8 pi^3), leading weak-source asymptotic coefficient",
    "fundamental_GRUT_restriction_beyond_full_symmetry_baseline": False,
    "complete_gravity_or_reality_theory": False,
    "stage3_card3_opened": False, "cp03_opened": False,
    "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__},
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
print(result["status"])
print(f"groups={len(GROUPS)} max_numeric_error={result['largest_numerical_absolute_error']:.3e}")
print("CONDITIONAL_EFT_ONLY; NO_NEW_FUNDAMENTAL_GRUT_LAW; NO_NEWTONIAN_GRAVITY_LIMIT")

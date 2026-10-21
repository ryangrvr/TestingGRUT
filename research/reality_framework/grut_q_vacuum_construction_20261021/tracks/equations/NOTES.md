# Track `equations`: blind derivation from the brief's action

21 October 2026. Construction pass for issue #13 (BRIEF.md in the parent directory).

**Read level and blindness.** This derivation does not use the Klinkhamer-Volovik (KV) papers. None was opened. The
`literature` track extracts them, and this track is the blind cross-check. Everything here is derived from the action
the brief transcribes. Standard results from other literature (Brans-Dicke/f(R) equivalence, 4-form flux, fifth-force
bounds) are labelled [R] (recalled), with identifiers. No downloads.

**Labels.** [D] derived here; [C] computed here by code in `code/`, with output in `results/`; [R] recalled; [E] estimate;
[A] argued. Conventions follow the brief: signature (-,+,+,+), hbar = c = 1, K = M^2 = 1/(8 pi G),
M = 2.435e18 GeV.

## 0. Verdicts in one screen

| # | Result | Label |
|---|--------|-------|
| 1 | The brief's (q), (g) and rho_V equations are correct: signs, factors and the -mu q term. | [D][C] |
| 1 | Off-shell identity nabla^m E_mn = -q d_n mu(x), so matter is conserved iff mu is constant. It holds exactly on random non-FRW and FRW jets, and fails for a wrong sign or a dropped term. | [D][C] |
| 2 | For fixed mu, the theory is Brans-Dicke with omega_BD = 0 and phi = K(q), U_mu(phi) = eps(q) - mu q. It equals f(R) where invertible. mu labels a one-parameter family of such theories. | [D] |
| 3 | FRW first-order system; the constraint propagates as Cdot = -4 H C exactly. | [D][C] |
| 4 | m_q^2 = (2/3) K(q0) eps''(q0) / K'(q0)^2. No ghost iff K > 0; no tachyon iff eps'' > 0. G_eff(r) = G_*(1 + e^{-m r}/3), gamma = (1 - e^{-m r}/3)/(1 + e^{-m r}/3). | [D][C] |
| 5a | With mu = mu0[eps + c]: rho_V,eq = 0 for every c, dq0/dc = 1/(q0 eps''), dmu0/dc = 1/q0. An equilibrium exists only for c in the range of q eps' - eps on the stable branch. | [D] |
| 5b | With fixed numerical mu: rho_V,eq = dc + (4/3) dc^2/(K m_q^2) + O(dc^3), so c gravitates one-for-one. H^2 = dc/(3K0) + O(dc^3): the second order cancels in H. | [D][C] |
| 5c | mu is conserved through a transition. With pre-transition mu, rho_V,inf = -Delta + (4/3) Delta^2/(K m_q^2) + O(Delta^3): AdS-type, and flat FRW recollapses. With post-transition (= asymptotic) mu, rho_V,inf = 0 to all orders. A second-order term never survives alone as a constant; it appears only as the decaying, dust-like oscillation energy (4/3) Delta^2/(K m_q^2) (x/sinh x)^2. | [D][C] |
| 6 | Zero needs eps = the asymptotic, fully renormalized total. The selection rule is a condition on the conserved integration constant, i.e. on the state, and it refers to the final state. The dynamics does not enforce it. For m_q >> H the theory is GR plus an integration-constant Lambda (Henneaux-Teitelboim-like) plus O(H^2/m_q^2). | [D][A] |

## 1. Field equations

### 1.1 Variables [D]

S = ∫ d^4x sqrt(-g) [ (1/2) K(q) R - eps(q) ] + S_m[psi, g].
F = dA, F_{κλμν} = 4 ∂_[κ A_λμν]. Define F~ ≡ (1/24) [κλμν] F_{κλμν} = (1/6) [κλμν] ∂_κ A_λμν. F~ is a density and
does not depend on the metric. Then F_{κλμν} = q sqrt(-g) [κλμν] with q = F~/sqrt(-g), which is equivalent to
q^2 = -(1/24) F_{κλμν} F^{κλμν}.
Define ℒ(q, R) = (1/2) K(q) R - eps(q).

### 1.2 Variation in A: the q-equation and mu [D][C]

δ_A (sqrt(-g) ℒ) = ∂ℒ/∂q δF~ = (1/6) ∂ℒ/∂q [κλμν] ∂_κ δA_λμν. Integrating by parts gives
[κλμν] ∂_κ (∂ℒ/∂q) = 0 for all λμν, so ∂_κ(∂ℒ/∂q) = 0. Hence

    (q)   mu ≡ -∂ℒ/∂q = eps'(q) - (1/2) K'(q) R = constant over spacetime.

mu is an integration constant, not a Lagrangian parameter. It is fixed by the state or boundary data: it is the
"electric flux" of the 3-form. [R] Brown-Teitelboim (Nucl. Phys. B297, 787 (1988)) and Duncan-Jensen
(Nucl. Phys. B336, 100 (1990)) describe holding mu fixed at the boundary by adding -∫ mu F~; not reread.
Minisuperspace check (`code/frw_minisuperspace.py`) [C]: with ds^2 = -N^2dt^2 + a^2dx^2 and A_123 = A(t), the
momentum conjugate to A is exactly K'(q)R/2 - eps'(q), and it is conserved.

### 1.3 Variation in g [D][C]

At fixed A, δq = -(q/2) g^{μν} δg_μν. Collecting the δ sqrt(-g), δq and K δR pieces gives

    K G_μν = T_μν + [q eps' - eps - (1/2) q K' R] g_μν - (g_μν □ - ∇_μ∇_ν) K.

Using (q):

    (g)   K(q) G_μν = T^m_μν - (eps - mu q) g_μν - (g_μν □ - ∇_μ ∇_ν) K(q),
          rho_V ≡ eps(q) - mu q          (flat equilibrium, R = 0: eps' = mu, so rho_V = eps - q eps').

This agrees with the brief. Checks [C]:
- The minisuperspace N-equation equals minus a^3 times the Friedmann constraint below.
- The a-equation equals -3a^2 times the covariant space-space equation.
- The covariant (g) equation, reduced by an explicit Christoffel computation, gives the same two FRW equations with
  zero residual, for generic K, eps, a(t) and q(t).

All results are in `results/frw_minisuperspace.json`.

### 1.4 Noether / Bianchi identity [D][C]

For any metric and any q(x), define mu(x) := eps'(q) - K'(q)R/2 and
E_μν := K G_μν + (eps - mu(x) q) g_μν + (g_μν□ - ∇_μ∇_ν)K. Then

    ∇^μ E_μν = -q ∂_ν mu(x)                                   (identity)

Proof: ∇^μ(K G_μν) = G_νμ ∇^μ K, and ∂_ν(eps - mu q) = (1/2) R ∂_ν K - q ∂_ν mu. Also
∇^μ[(g□ - ∇∇)K]_μν = ∇_ν□K - □∇_νK = -R_νμ ∇^μ K. These terms sum to -q ∂_ν mu.

On shell, E = T^m, so ∇^μ T^m_μν = -q ∂_ν mu. The q-equation (∂mu = 0) is therefore equivalent to matter
conservation where q ≠ 0. Equivalently, if E is built with a constant mu, ∇^μ E_μν = (mu(x) - mu) ∂_ν q.

Checks:
- [C] Symbolic FRW with generic K, eps, a(t), q(t): residual 0 (`frw_minisuperspace.py`).
- [C] Exact rational-arithmetic jets (`code/bianchi_jet_check.py`): metrics with random non-diagonal entries depending
  on all four coordinates (3 configurations), and random FRW (2), with random cubic K(q), eps(q) and q(x). The
  identity holds exactly (residual 0). It fails, as it should, when the sign of mu q is flipped or the K-derivative
  terms are dropped. The constant-mu form holds exactly. See `results/bianchi_jet_check.json`.

### 1.5 Dual action and degrees of freedom [D]

For constant mu, ∫ mu F~ d^4x = ∫ sqrt(-g) mu q is a total derivative. So S and
S_dual = ∫ sqrt(-g)[ K(q) R/2 - eps(q) + mu q ] + S_m, with q varied as an independent scalar, give identical field
equations: δq reproduces (q), and δg reproduces (g) with no δq(g) term.

Degrees of freedom:
- 2 graviton polarizations.
- 1 scalar mode (the "q-mode") if K'(q0) ≠ 0.
- 0 local degrees of freedom from A, since a 3-form in 4D does not propagate.
- One global constant mu per connected spacetime. It changes only through membrane nucleation, which is outside this
  action ([R] Brown-Teitelboim 1987/88).

If K' ≡ 0, q is locked algebraically by eps'(q) = mu, and the theory is exactly GR with
Lambda_eff = eps(q_mu) - mu q_mu = -eps*(mu), the Legendre transform. This is the 4-form-flux / Henneaux-Teitelboim
picture of the cosmological constant as an integration constant ([R] Aurilia-Nicolai-Townsend 1980;
Duff-van Nieuwenhuizen 1980; Hawking 1984; Henneaux-Teitelboim, Phys. Lett. B222, 195 (1989)).

## 2. Trace equation, q dynamics, scalar-tensor map

### 2.1 Trace and the second-order equation [D][C]

The trace of (g) is -K R = T - 4(eps - mu q) - 3□K. Eliminating R with (q), R = 2(eps' - mu)/K', gives

    3 □K(q) = T - 4(eps - mu q) + 2 K (eps' - mu)/K'
    3 K' □q + 3 K'' (∇q)^2 = T - 4(eps - mu q) + 2 K (eps' - mu)/K'

Equivalently □phi = (2/3)(phi U' - 2U) + T/3 with phi = K and U as in 2.2 (checked in
`code/linear_and_energy_origin.py`).

### 2.2 The map [D]

On any interval where K'(q) ≠ 0, set phi = K(q) and q = q(phi). For a given mu the field equations are exactly those
of

    S_BD = ∫ sqrt(-g) [ (1/2) phi R - U_mu(phi) ] + S_m,   U_mu(phi) = eps(q(phi)) - mu q(phi),   omega_BD = 0,

with U' = (eps' - mu)/K' and U'' = [K' eps'' - (eps' - mu) K'']/K'^3. In the usual Brans-Dicke normalization
(16π)^-1 ∫ sqrt(-g)[Φ R - (ω/Φ)(∂Φ)^2 - V], this is Φ = 8π K(q), V = 16π U_mu, ω = 0. In this map:
- The 3-form equation becomes the phi-equation R = 2U'.
- mu becomes a parameter of the potential.

So q-theory is a one-parameter family of omega = 0 BD theories, and the state chooses the member.

Where in addition eps'' - K''R/2 ≠ 0, eliminating q gives an f(R) theory, S = ∫ sqrt(-g) f(R)/2, with
f(R) = K(q_R)R - 2eps(q_R) + 2 mu q_R and q_R defined by (q). Then:
- f'(R) = K(q_R), by the envelope theorem.
- f''(R) = K'^2 / [2(eps'' - K''R/2)].
- The standard scalaron mass (f' - R f'')/(3f'') at R = 0 reproduces m_q^2 of section 4.

[R] The omega = 0 BD ≡ f(R) equivalence is due to O'Hanlon (1972) and Teyssandier-Tourrenc (1983); for a review see
Sotiriou-Faraoni, Rev. Mod. Phys. 82, 451 (2010), arXiv:0805.1726. If KV state this map, the credit is theirs; this
track did not check whether they do.

## 3. Flat FRW with dust, radiation and a matter-sector vacuum term

### 3.1 Equations [D][C]

ds^2 = -dt^2 + a^2 dx^2, q = q(t), K = K(q), Kdot = K' qdot. Matter is
rho_m = rho_d + rho_r + V_m and p_m = rho_r/3 - V_m.

    Friedmann (00):   3 H^2 K + 3 H Kdot = rho_m + (eps - mu q)
    space-space:      -(2 Hdot + 3H^2) K = p_m - (eps - mu q) + Kddot + 2 H Kdot
    q-equation:       eps'(q) - 3 K'(q)(Hdot + 2H^2) = mu
    matter:           rho_d' = -3H rho_d ;  rho_r' = -4H rho_r - dV_m/dt  (latent heat goes into radiation)

V_m enters only T^m and does not depend on q. In the metric equation it acts exactly like eps -> eps + V_m, and it is
absent from (q) [D].

Supplied history [P]: V_m(N), N = ln a, a tanh step of depth Delta and width w. A temperature-based version is
possible: with rho_r = a_r T^4, Tdot = -4H rho_r/(4 a_r T^3 + dV_m/dT). This conserves energy but is a bag-model
simplification; the thermodynamically complete form uses p(T) = a_r T^4/3 - V_m(T) and rho = T p' - p [A].

### 3.2 First-order system for numerics [D]

Use y = (q, v = qdot, H, rho_d, rho_r [, N]) and R_q(q) ≡ 2(eps'(q) - mu)/K'(q):

    qdot = v
    vdot = { -K'' v^2 - 3H K' v + [rho_d + 4 V_m + 4(eps - mu q) - K R_q]/3 } / K'
    Hdot = R_q/6 - 2H^2
    rho_d' = -3H rho_d,  rho_r' = -4H rho_r - V_m',  N' = H

Initial data must satisfy C ≡ 3H^2K + 3HK'v - rho_m - (eps - mu q) = 0 (solve for H). For q0 take the quasi-static
root of the trace equation; the initial v is the quasi-static drift.

### 3.3 Constraint propagation [C]

Symbolically, for generic K, eps, arbitrary V_m(t) and the system above, Cdot = -4 H C exactly. So the constraint
surface is invariant: a numerical violation decays as a^-4 during expansion and grows during collapse. In the toy runs
the maximum relative violation is at most 1.3e-9.

### 3.4 Quasi-static reduction for the physical hierarchy [D]

For m_q >> H, q = q_qs + (oscillation) + O(H^2/m_q^2), with K' δq_qs = K'^2 (rho_d + 4 dc)/(2 K eps'') and
dc ≡ V_m(t) - c*(mu). Here c*(mu) is the matter vacuum value for which mu is the Minkowski equilibrium value. The
Friedmann equation becomes

    3 K0 H^2 = rho_d + rho_r + dc(t) + O(dc^2/(K m_q^2), rho_d^2/(K m_q^2))

This is GR with a vacuum energy V_m(t) - c*(mu): an integration-constant Lambda, as in section 1.5. The q-mode shows
up only at O(H^2/m_q^2), plus oscillations excited at transitions (section 5.4).

## 4. Linearization about Minkowski (q0, mu0)

### 4.1 Mass [D][C]

Perturb at fixed mu = mu0, since a change of mu is a change of state and not a fluctuation:

    (□ - m_q^2) δK = T/3,     m_q^2 = (2/3) K(q0) eps''(q0) / K'(q0)^2 = (2/3) K(q0) / [ χ (q0 K'(q0))^2 ],

with χ^-1 ≡ q0^2 eps''(q0), my own "vacuum compressibility" normalization; I have not checked KV's. The mass depends
on q0 only through K, eps'' and K' evaluated there. Checks:
- The Einstein-frame mass (2/3) phi0 U''(phi0) agrees, with sigma = sqrt(3/2) M ln(phi/M^2) and
  V = M^4 U/phi^2 (`linear_and_energy_origin.json`).
- The toy FRW oscillation frequency is 4.0824828, against the predicted 4.0824829 [C].

The mode disappears (m_q → ∞) as K' → 0.

### 4.2 No ghost / no tachyon [D]

The Einstein-frame kinetic coefficient is ∝ (omega_BD + 3/2) = 3/2 > 0, and the graviton is healthy if K(q0) > 0.
So there is no ghost iff K(q0) > 0, and no tachyon iff eps''(q0) > 0 (given K > 0). A propagating mode also requires
K'(q0) ≠ 0.

De Sitter states with H << m_q are stable [A]. The standard f(R) criterion is (f' - R f'')/f'' > 0 at R_dS ([R]
Sotiriou-Faraoni). The toy dS run (5.4, case C) converges [C].

### 4.3 Newtonian limit, Yukawa, PPN [D][C]

Take ds^2 = -(1+2Ψ)dt^2 + (1-2Φ)dx^2. The linearized equations are:
- 2 K0 ∇^2Φ = rho + ∇^2 δK
- Ψ - Φ = -δK/K0
- (∇^2 - m^2) δK = -rho/3

The linearized Einstein tensor components were checked symbolically. For a point mass M_p:

    δK = M_p e^{-m r}/(12π r),  Ψ = -(G_* M_p/r)(1 + e^{-m r}/3),  Φ = -(G_* M_p/r)(1 - e^{-m r}/3),  G_* = 1/(8π K(q0))
    gamma_PPN(r) = (1 - e^{-m r}/3)/(1 + e^{-m r}/3);   Yukawa range lambda = 1/m_q, strength alpha = 1/3
    quasi-static Fourier form: G_eff(k) = G_* [1 + k^2/(3(k^2 + m_q^2))]

The ODEs and flux normalizations were verified symbolically.
- Light q (m r << 1): gamma = 1/2 and G_lab = (4/3) G_*. Cassini requires |gamma - 1| < 2.3e-5 ([R] Bertotti et al.,
  Nature 425, 374 (2003)), which excludes this unless m_q R_sun ≳ 10, i.e. m_q ≳ 3e-15 eV. Short-range torsion
  balances (alpha = 1/3) need lambda ≲ few × 10 μm, i.e. m_q ≳ 4 meV ([R], order of magnitude; Eöt-Wash, Lee et al.,
  PRL 124, 101101 (2020)) [C].
- Heavy q: G_lab = G_*, gamma = 1.

The cosmological G in the Friedmann equation is G_cos = 1/(8π K(q_bg)) = G_*, up to δK/K = rho_m/(3 K m_q^2)
≈ H^2 Ω_m/m_q^2. That is 1.6e-121 today for natural m_q [C]. So Cavendish G equals cosmological G unless q is light.

## 5. Energy-origin test: eps -> eps + c

### 5.1 (a) State given by mu = mu0[eps + c] [D][C]

The equilibrium conditions are eps'(q0) = mu0 and eps(q0) + c - mu0 q0 = 0. So rho_V,eq = 0 for every c, as an exact
statement at any order in c. Implicit differentiation gives:

    dq0/dc = 1/(q0 eps''(q0)),        dmu0/dc = 1/q0,
    d2q0/dc2 = -(eps'' + q0 eps''')/(q0 eps'')^3,   d2mu0/dc2 = -1/(q0^3 eps''),
    dK(q0)/dc = K'(q0)/(q0 eps''(q0)) = 2K/(3 m_q^2 q0 K')   (G_eq shifts; with natural scales δK/K ~ c/M^4).

The constant is absorbed by moving q0. Toy check: the finite difference of mu0 equals the secant of 1/q0 [C].

Existence [D]: the equilibrium needs c = P(q0) with P(q) ≡ q eps' - eps, and P' = q eps''. On a stable (eps'' > 0)
branch, P is monotonic in |q|. For a convex branch through q = 0 this requires eps(0) + c > 0. A negative total
vacuum energy at q = 0 has no stable Minkowski equilibrium on such a branch. The mechanism compensates one sign of
total vacuum energy unless eps is non-convex.

### 5.2 (b) Fixed numerical mu [D][C]

Static (dS/AdS) states satisfy 3H^2K = rho_V, R = 4rho_V/K, eps'(q) - mu = 2K'(q) rho_V/K(q), and
rho_V = eps + c - mu q. Define c* by mu = mu0[eps + c*] and dc = c - c*. Expanding with generic Taylor data
(K0, K1, K2, K3, e2, e3, e4) at q0:

    q_s - q0 = 2K1 dc/(K0 e2) + O(dc^2)
    rho_V,eq = dc + [2 K1^2/(K0^2 e2)] dc^2 + [8 K1^2 (3 K2 e2 - K1 e3)/(3 K0^3 e2^3)] dc^3 + O(dc^4)
             = dc + (4/3) dc^2/(K0 m_q^2) + ...
    H_s^2   = dc/(3 K0) + [2 K1^2 (3 K2 e2 - K1 e3)/(9 K0^4 e2^3)] dc^3 + O(dc^4)   (no dc^2 term)

So d rho_V,eq/dc = 1 at the Minkowski point exactly: a fixed mu lets any constant gravitate one-for-one. With
K' ≡ 0, rho_V,eq = c - eps*(mu) exactly. With the local G = 1/(8π K(q_s)), the observable is rho_V. Measured against
K0, the expansion rate sees only dc, through second order.

### 5.3 (c) A transition: matter vacuum drops by Delta, mu conserved [D]

mu is conserved through the transition: (q) holds everywhere and V_m does not involve q. If V_m(q) did depend on q,
∂V_m/∂q would join the definition of mu, which would still be conserved. Because q obeys a second-order equation,
q and R stay continuous even through a sudden drop, unlike GR, where R = -T/K jumps.

The late-time state is a fixed point of 3.2 with matter diluted: Minkowski iff mu = mu0[eps + c_inf]. Otherwise it is
the dS/AdS root of 5.2 with dc = c_inf - c*(mu). That gives:

| choice of mu | dc | rho_V,inf | attractor |
|---|---|---|---|
| mu0 of the pre-transition total eps + c_before | -Delta | -Delta + (4/3)Delta^2/(K m_q^2) - c3 Delta^3 | AdS root; flat FRW turns around and crunches |
| mu0 of the post-transition total (= asymptotic if this is the last transition) | 0 | 0 exactly, all orders | Minkowski |
| mu0 of the asymptotic total after later transitions Delta_i | 0 | 0 exactly | Minkowski |
| mu0 of an earlier equilibrium, missing later Delta_i | -Σ Delta_i | -ΣΔ + (4/3)(ΣΔ)^2/(K m_q^2) | AdS / crunch |
| mu0 of the bare eps (matter vacuum excluded) | c_inf | c_inf + (4/3) c_inf^2/(K m_q^2) | dS or AdS |
| vacuum energy rises by Delta (pre mu) | +Delta | Delta + (4/3)Delta^2/(K m_q^2) | de Sitter |

Here c3 = 8K1^2(3K2e2 - K1e3)/(3K0^3e2^3).

**Which order survives [D]:** the first order survives unless mu refers to the asymptotic total eps. When it does,
every order vanishes. Within this action, a second-order term never survives on its own as a constant remnant. Under
post-transition or asymptotic mu it appears only as transients:
1. The dust-like oscillation energy of the q-mode after the transition,
   rho_osc = (4/3) Delta^2/(K m_q^2) × (x/sinh x)^2 with x = π m_q τ/2, where τ is the transition duration
   [D: linear oscillator driven by a tanh step]. It redshifts as a^-3.
2. The pre-transition dS offset, rho_V,pre = Delta + (4/3)Delta^2/(K m_q^2).
3. A matter-induced rho_V ≈ rho_d^2/(12 K m_q^2) → 0.

### 5.4 Toy numerical confirmation [C]

`code/frw_transition_toy.py`, `results/frw_transition_toy.json`. Toy inputs [P]: eps = 1/2 + q^2/2,
K = 1 + 0.2(q-1), so m_q^2 = 16.67; Delta = 0.02, rho_r0 = 0.3, rho_d0 = 0.01; transition at ln a = 0.8; widths
w = 0.003, 0.1, 0.5. Units K(q0) = 1. In this toy m_q/H ~ 40 at the transition; physically it is ~1e32 (section 7).

- **A (post-transition mu).** Minkowski attractor. rho_V,tot falls to ≤1e-8 by t = 3000 (H = 1.8e-4).
  - Sudden transition (w = 0.003): oscillation energy × a^3 = 3.493e-4, against 3.488e-4 predicted
    ((4/3)Δ^2/(K m^2) a_c^3 × form factor 0.989), a 0.1% match.
  - w = 0.1: suppression 7.1e-4, against a form factor of 5.4e-4. The ~30% gap comes from H varying during the
    transition.
  - w = 0.5: suppression ≤1e-5. The floor is set by the background, not by the form factor.
  - Control with Delta = 0: 5e-10.
  - The q oscillation frequency matches m_q to 4e-7.
- **B (pre-transition mu, drop).** The static AdS root is rho_V = -0.0199683, matching the series
  -Δ + (4/3)Δ^2/(Km^2) to 1e-16, since c3 = 0 for this toy. Dynamically, rho_V = -0.019973 at turnaround (t ≈ 10.7),
  followed by recollapse.
- **C (vacuum rises by Δ, pre mu).** de Sitter attractor: rho_V,late = 0.020032 = the static root
  = Δ + (4/3)Δ^2/(K m^2); H^2 = 0.0066667 = Δ/(3K0), with no second-order term in H^2. Confirmed.

## 6. What this implies for "which eps does the equilibrium condition refer to" [D][A]

1. **The zero requires the asymptotic, fully renormalized total.** rho_V,inf = 0 iff mu = mu0[eps_tot,inf], where
   eps_tot,inf(q) = eps(q) + (every q-independent vacuum energy of the matter sector in the final state: all
   thresholds, loops and transition remnants, evaluated at the true final minimum). Any other choice leaves the
   missing contributions at first order (5.3 table). The bare eps, pre-transition eps, or any eps that omits later
   transitions are all excluded as readings of "the equilibrium of the whole" if the residual is to vanish [D].
2. **The condition is on the state, not on the dynamics.** mu is a conserved Gauss-law constant (1.2). Transitions
   in sectors that carry no 3-form charge cannot change it, and the q-dynamics (K' ≠ 0) relaxes q towards the
   equilibrium set by mu, not mu towards zero vacuum energy. Imposed at early times ("the universe began in
   equilibrium"), the condition gives -ΣΔ_i; for the electroweak crossover alone that is -4.7e54 ρ_Λ,obs. To get zero,
   mu must be chosen with knowledge of the final vacuum: a global or final-state boundary condition [D/A].
3. **It passes the energy-origin test by construction, as a selection rule.** Under eps -> eps + c the rule
   re-selects mu (dmu0/dc = 1/q0) and keeps rho_V = 0 (5.1). A fixed number mu does not (5.2). No Lagrangian parameter
   is retuned, but the selected number mu moves by δc/q0 for every heavy-threshold contribution δc. The tuning is
   carried by the state, as in unimodular or Henneaux-Teitelboim gravity (1.5, 3.4) [D].
4. **It fixes the residual to exactly zero, not to ρ_Λ.** Matching ρ_Λ,obs needs mu = mu0[eps_tot,inf] - ρ_Λ/q0 + O(ρ_Λ^2):
   a calibrated offset of relative size ρ_Λ/M^4 = 7.2e-121 when q0 mu0 ~ M^4 [D/C]. The remaining equilibrium-branch
   sources are tiny or excluded:
   - The matter-induced ρ_V ≈ ρ_m^2/(12 K m_q^2) is 5e-169 GeV^4 for natural m_q.
   - Raising it to ρ_Λ needs m_q ≈ 0.2 H0. That is outside the regime m_q >> H of the formula, and its gamma = 1/2 is
     excluded locally unless screened [C/E].
   - The second-order (4/3)Δ_EW^2/(K m_q^2) = 8e-58 GeV^4 = 3e-11 ρ_Λ is not a constant remnant in any reading (5.3).
5. **Proposed GRUT formalization, made precise.** "The closed whole has no exterior" ⇒ "the state of the whole has
   mu = mu0[eps_tot,inf]", i.e. vanishing vacuum pressure of the final, asymptotic total. It is an assumption on the
   state, credited to KV for the equilibrium idea [A]. Whether KV also make it refer to the final state is for the
   literature track to settle. As derived here, it removes every constant (sign and scale: zero). It is c-covariant.
   Its consequences are those of GR with zero final Λ plus a heavy scalar, so it predicts no held-out quantity at the
   current level.

## 7. Realistic-scale numbers [C/E]

From `code/realistic_numbers.py` and `results/realistic_numbers.json`. Inputs: M and ρ_Λ from the brief; m_H = 125.25,
v = 246.22 GeV [R]; bag constant B^1/4 = 0.145-0.235 GeV [R]; H0 = 67.4, Ω_m = 0.315 [R]; natural q-scales [E].

| Quantity | Value |
|---|---|
| natural m_q (q0^2 eps'' = M^4, q0 K' = M^2) | sqrt(2/3) M = 2.0e18 GeV [E] |
| Δ_EW = m_H^2 v^2/8 | 1.19e8 GeV^4 (tree level) [D/C] |
| Δ_QCD | 4.4e-4 to 3.0e-3 GeV^4 [R] |
| Pre-transition-mu remnant | -4.7e54 ρ_Λ (EW); -1.8e43 to -1.2e44 ρ_Λ (QCD) |
| (4/3)Δ_EW^2/(K m_q^2) | 8.0e-58 GeV^4 = 3.2e-11 ρ_Λ. It would equal ρ_Λ only at m_q = 1.1e13 GeV, and even then it is not a constant remnant. |
| m_q/H at T = 100 GeV | 1.4e32, so x/sinh x is effectively 0: no oscillation energy is deposited |
| matter-induced δG/G today | 1.6e-121 |

## 8. Inputs and provenance (this track)

- Action and conventions: from the brief, [P] supplied, borrowed from KV and credited to them.
- K(q), eps(q): unspecified generic functions; results are stated through their Taylor data at q0.
- Toy functions and numbers in 5.4: supplied [P], illustrative, not calibrated.
- mu: integration constant (the state). Its selection is the open assumption (section 6).
- No parameter was calibrated to ρ_Λ,obs. Section 6.4 states what a calibration would fix.

## 9. Limits and open items

- Classical field equations only; no quantum stability of the q-sector or loop corrections to K(q) and eps(q) beyond
  the constant c.
- The adiabatic form factor is exact only for a constant-frequency oscillator. It is within 30% at w = 0.1 and fails
  as a predictor at w = 0.5 (H-variation floor).
- Inhomogeneous relaxation (bubbles, domain walls at a first-order transition) was not studied.
- Membrane nucleation (Brown-Teitelboim), which changes mu, lies outside the action.
- Cross-check against the literature track still needed:
  - KV's form of the q-equation and of m_q, including their compressibility normalization;
  - how KV fix mu: pre-transition, post-transition or asymptotic;
  - whether KV's "kicked vacuum" Δ^2/M^4-type remnant (recalled scaling only, not read) arises from a different
    setup. In this action it does not survive as a constant.

## Files

- `code/frw_minisuperspace.py` → `results/frw_minisuperspace.json`: sections 1.2, 1.3, 1.4 (FRW), 3.3.
- `code/bianchi_jet_check.py` → `results/bianchi_jet_check.json`: section 1.4.
- `code/linear_and_energy_origin.py` → `results/linear_and_energy_origin.json`: sections 2, 4, 5.1, 5.2, 5.3.
- `code/frw_transition_toy.py` → `results/frw_transition_toy.json`: section 5.4.
- `code/realistic_numbers.py` → `results/realistic_numbers.json`: section 7.

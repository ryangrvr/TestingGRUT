# Track `obstruction`: attempt to refute the energy-origin obstruction, and its exact scope

21 October 2026. Brief: `../../BRIEF.md`. Object of study: Codex's energy-origin obstruction
(`git show 02b4515:research/reality_framework/grut_scale_feasibility_20261021/NEXT_DIRECTION.md`).
Labels: **[D]** derived here, **[C]** computed by `code/obstruction_checks.py` (record:
`results/obstruction_checks.json`, 34/34 checks pass), **[Q]** quoted from a source at the read level
given in §9, **[E]** order-of-magnitude estimate, **[A]** argued (not a theorem). Nothing here is a
new theory, an empirical validation or a confirmatory test. Numbers built from H0 and Omega_Lambda
are **retrospective** (already-exposed data). Conventions: signature (-,+,+,+); reduced Planck mass
M = (8 pi G)^(-1/2) = 2.435e18 GeV; hbar = c = 1.

## 0. Result

1. **The obstruction is correct as scoped [D].** Every closure property in FORMAL_CORE §§1-5
   (global law, restriction, reduced dynamics, NZ/TCL, influence functional, Gibbs complements,
   Sigma = I + D, KMS, modular flow, athermality budget, fluctuation relations, internal-clock
   relations) is invariant under H -> H + cI (§1, checks A1-A23). The only non-invariant items are
   labels (<H>, Z, F) and rules that state an absolute energy as a supplied number (a constraint
   value, a microcanonical shell at a stated E, a renormalization condition).
2. **Where the origin does become physical, it is supplied from outside closure [D/A].** It is
   physical whenever a dynamical variable couples to the constant: an internal control whose energy
   offset depends on its coordinate feels a force (Casimir structure, A19-A20); a pointer coupled to
   H as an operator reads c (A23); the metric couples through sqrt(-g), so a vacuum density is an
   energy that depends on a dynamical coordinate (§1.5, D5). The origin is *fixed* only by a symmetry
   or boundary statement: the Poincare algebra fixes the global generator (H|0> = 0, B1-B2) but not
   the c-number in <T_mu nu>; time reparametrization fixes the Hamiltonian-constraint value at zero
   but leaves Lambda a parameter that the constraint relates to the state (§2, C1-C3, C7-C8);
   asymptotic flatness fixes Lambda_tot = 0 by boundary condition; Wald's Minkowski axiom fixes the
   matter side of the unphysical split. None selects the observed residual.
3. **Two precision corrections [D/Q].** (a) The task text says Wald's axioms leave c1 g_mu_nu. That
   holds for the Hollands-Wald locally covariant scheme (c-number counterterm a m^4) and for Wald's
   axioms without the Minkowski-vacuum normalization. With that axiom, the matter-side c1 is fixed by
   convention. This fixes the split, not the total, because the gravitational Lambda_b stays free
   (§3). (b) Unimodular gravity does not refute the obstruction. It moves the free constant from a
   coupling to a conserved quantity of the global state, the conjugate of 4-volume (§4). The
   obstruction then becomes a state-selection problem, the slot where earlier GRUT closure work found
   no rule.
4. **Minimal escape class [A]** (§6.2). A principle that fixes the total renormalized rho_vac
   without a bare value statement needs at least one of these: (G) a global or nonlocal functional
   sensitive to the constant volume mode; (B) a background or vacuum-symmetry-breaking structure
   that removes that mode from the local metric equation; (S) a state or measure over a space of
   Lambda values; (X) degrees of freedom beyond a finite local 4D EFT. To also fix the sign and
   scale and stay stable, it further needs (R) a coupling-independent scale and (P) protection of
   the constant mode at all orders. Each item negates one stated hypothesis. Every class has
   published representatives; none is GRUT-native.
5. **GRUT closure supplies none of G, B, S, X, R, P by itself [A]** (§6.3).

## 1. Fixed-background quantum mechanics

### 1.1 The displayed identities [D; C: A1-A9]

With [H, cI] = 0: e^{-i(H+c)t} = e^{-ict} e^{-iHt}. Hence U_c rho U_c^dag = U rho U^dag,
[H_c, rho] = [H, rho], Tr_{Abar} of the evolved state is unchanged, Born probabilities Tr(rho E) are
unchanged, and spec(H_c) = spec(H) + c, so every gap E_k - E_l is unchanged. For the Gibbs state:
e^{-beta H_c} / Tr e^{-beta H_c} = e^{-beta c} e^{-beta H} / (e^{-beta c} Z). Thus Z_c = e^{-beta c} Z and
F_c = F + c. The partition function and free energy are labels that change; the state does not.
The code checks a 2x3 composite with a random interacting H, c = 7.31, t = 2.7, beta = 0.8:
deviations are <= 6e-15 (tolerance 1e-11), and F_c - F = c to 9e-16.

In-in versus in-out [D; C: A6]: <psi|U_c|psi> = e^{-ict} <psi|U|psi>. The vacuum energy is a phase of
the in-out amplitude. In-in (closed-time-path) objects contain U_c ... U_c^dag, so the phases cancel.
This is why the influence functional and S_CTP of FORMAL_CORE THM 2.5 are blind to c.

### 1.2 Every FORMAL_CORE closure property is shift-invariant [D]

| FORMAL_CORE item | how H enters | invariant under H -> H + cI? |
|---|---|---|
| DEF 1.1 law; THM 2.1 spectrum of rho, S, purity | Ad U(t) | yes (A2). <f(H)> is conserved for every f; the value <H> is a label |
| DEF 1.4 restriction; THM 2.6-2.7 marginals | through rho(t) | yes (A3) |
| THM 2.2 reduced equation | commutators; cI = c_A I x I + I x c_Abar I | yes (A4) |
| THM 2.3 NZ/TCL; 2.4 maps, divisibility | L = -i[H, .] | yes |
| THM 2.5 influence functional, S_CTP | forward and backward branches | yes: phases cancel (A6) |
| DEF 1.5 / THM 2.8 records, process tensor | internal unitaries; external controls | yes for states. A control-dependent c-number c(lambda) shifts work labels only (A18) |
| THM 2.9 entropies, recurrence | rho(t) | yes |
| COND 2.10 equilibration | gaps, d_eff | yes (A5) |
| COND 2.11 typicality | measure on a subspace H_R | yes if H_R is named by a spectral projector; **no** if it is named by an absolute energy value. That value is then supplied |
| COND 2.12 Lieb-Robinson; 2.13 WAY | commutators, conservation | yes |
| THM 2.14 athermality budget; 2.15 protect-or-form | stationary sigma, [H, O] | yes (A7) |
| COND 2.16 Sigma = Delta S_A + beta Delta E_Abar = I + D | Delta E_Abar = Tr H_Abar (rho_Abar(t) - gamma) | yes, since Tr(rho_Abar(t) - gamma) = 0 (A14-A15) |
| COND 2.17 KMS reservoirs (infinite) | derivation delta(A) = i[H, A] | yes. The GNS generator's origin is fixed by the convention L Omega = 0 (standard Liouvillean) |
| COND 2.18 FDT; 2.19 chaos bound | spectral functions, OTOCs | yes |
| COND 2.20 Pauli / Unruh-Wald; speed limits | H bounded below; <H> - E_0; Delta H | yes |
| COND 5.1-5.3 type III, modular flow, natural states | state only | yes (§1.3) |

### 1.3 KMS and modular Hamiltonians [D; C: A10-A13]

KMS at beta is defined relative to the automorphism group alpha_t(A) = e^{iHt} A e^{-iHt}. Since
alpha_t^{(c)} = alpha_t identically, the KMS states of H and H + cI coincide (A10-A11, residual
<= 7e-15). The modular operator Delta_omega and the flow sigma_t^omega depend only on the faithful
state. In the finite standard form, log Delta = log rho_A x 1 - 1 x log rho_A' has no free constant.
The one-sided K_A = -log rho_A is unique once Tr rho_A = 1 (A12). Only the parametrization
rho_A = e^{-K}/Z leaves K defined up to K + k, with k absorbed in Z (A13). For a Gibbs restriction,
K_A = beta(H_A - F_A): the additive constant of the modular Hamiltonian is beta times the free
energy, which is exactly where c sits. It is fixed by normalization and carries no dynamics. In
type III algebras (COND 5.1) there is no one-sided K_A, only log Delta and relative modular
operators. **Answer:** the modular flow, the only modular object that enters KMS or the dynamics,
is independent of any constant. No modular or KMS property can select c.

### 1.4 Entropy production and fluctuation relations [D; C: A14-A18]

Sigma(t) = I(A:Abar)(t) + D(rho_Abar(t) || gamma_Abar) holds to 3e-15 and is unchanged when
H_Abar -> H_Abar + cI. Two-point-measurement work is W = E_m^(1) - E_n^(0). Under a constant c, W
is unchanged. Under a control-dependent c-number c(lambda), W and Delta F both shift by
c(lambda_1) - c(lambda_0), and <e^{-beta W}> = e^{-beta Delta F} holds in every case (A16-A18; the
computed shifts equal 3 to 1e-15).

### 1.5 Lemma: a constant becomes physical only when it depends on a dynamical variable [D; C: A19-A20]

Let the control be a dynamical coordinate x inside the whole:
H_tot = H_S x 1 + 1 x H_X + 1 x c(x_hat). Then c(x_hat) is an operator, not a multiple of the
identity, and exerts the force -c'(x) on the control. Numerically, an x-independent c leaves
<x(t)> unchanged (deviation 6e-13). With c = kappa x, the harmonic control is displaced by
-kappa(1 - cos wt)/(M w^2) (computed -0.7428 against -0.7434; grid error 7e-4). This is the
structure of the Casimir effect: differences of vacuum energy between control settings are physical,
a setting-independent constant is not. **Gravity is the case in which the "control" is the local
volume element sqrt(-g) at every point:** S_m -> S_m - C ∫ sqrt(-g) gives
Delta T_mu_nu = -C g_mu_nu (D5, sign checked symbolically). A vacuum density is therefore an energy
that depends on a dynamical coordinate. This is the precise reason the fixed-background argument
cannot be extended to the gravitating term, as Codex states.

### 1.6 Internal clocks and time-energy relations [D; C: A21-A22]

Page-Wootters (DEF 1.6; a supplied construction): (H_S + H_C) Psi = eps Psi. Shifting
H_S -> H_S + c leaves the physical subspace unchanged **iff** eps -> eps + c. In the computed case
the subspace dimensions are 2, then 0 with eps held fixed, then 2 with eps shifted. Relational
states <tau|_C Psi> are unchanged under the combined shift. With an internal clock, the origin
reappears only as the constraint value eps. Nothing in the finite framework fixes eps; in GR,
reparametrization invariance does (§2). Mandelstam-Tamm (Delta H), Margolus-Levitin (<H> - E_0),
Pauli and Unruh-Wald (H bounded below) are all invariant.

### 1.7 Apparatus coupled to the absolute energy [D; C: A23]

The interaction g(H_S + c) x p_pointer differs from g H_S x p_pointer by g c (1 x p_pointer), which
is an operator on the pointer. The pointer moves by g c t (0.3 computed and predicted, to 1e-6). The
origin is set by the supplied interaction term, that is, by which operator the apparatus couples to.
Gravity's coupling to T_mu_nu has the same structure.

### 1.8 A symmetry algebra fixes the global generator, not the local density [D; C: B1-B2]

In 1+1, with H = i d_t, P = -i d_x, K = i(x d_t + t d_x): [K, P] = iH, and [K, P] = i(H + c)
fails unless c = 0. A Poincare-invariant vacuum (P|0> = K|0> = 0) therefore has
H|0> = -i[K, P]|0> = 0. The 3+1 Poincare group has no nontrivial central extension (Bargmann), so
the global generator's origin is fixed by symmetry. However, <0|T_mu_nu(x)|0> = -rho eta_mu_nu is
Poincare-invariant for every rho, and the generator is ∫(T_00 - <T_00>) d^3x. The algebra constrains
only the subtracted operator. The same holds for SO(4,1) and <T> = -rho g in de Sitter. **Gravity
couples to the local density, which the algebra leaves free [A].**

### 1.9 Answer to (1)

No closure property in the framework is non-invariant. The obstruction holds on its stated scope
[D]. What would break invariance is extra structure: a dynamical coupling to the constant, or a
principle that states an absolute value.

## 2. Closed gravitational systems

### 2.1 ADM form and how Lambda enters [D]

S = ∫ sqrt(-g) (M^2/2)(R - 2 Lambda_b) + S_m. In ADM variables,
H_grav = ∫_Sigma d^3x (N calH + N^i calH_i) (+ surface terms when Sigma is noncompact), with

  calH = (2/(M^2 sqrt h)) (pi^ij pi_ij - pi^2/2) - (M^2/2) sqrt h (^3R) + sqrt h (rho_m + rho_Lambda,b) ≈ 0,
  rho_Lambda,b ≡ M^2 Lambda_b.

The lapse is a Lagrange multiplier, so the constraint value is **zero**, fixed by time
reparametrization invariance. It is not a free eps as in §1.6.

* Relabelling (rho_m,vac + C, rho_Lambda,b - C) leaves calH unchanged (C1, symbolic).
* The uncompensated shift changes calH by C sqrt(h). In FRW the constraint changes by -a^3 C (C2).
  **A closed gravitational system does see the energy origin**, as Codex's scope note anticipated.
* What the constraint does with Lambda (C3): in FRW,
  rho_Lambda = 3M^2 (H^2 + k/a^2) - rho_m. This is how Lambda is *measured* (retrospectively) from
  state data. It is not predicted: solutions exist for every rho_Lambda, with H adjusted
  accordingly. **The constraint relates Lambda to the state; it does not fix Lambda [D].**
* Compact Sigma (spatially closed): there is no surface term, and H_total ≈ 0 on every solution for
  every Lambda_b. "The total energy of a closed universe is zero" imposes no condition on Lambda_b.
* Asymptotically flat: H = ∫(N calH + N^i calH_i) + E_ADM, with E_ADM >= 0 (Schoen-Yau; Witten) and
  Minkowski at zero. Asymptotic flatness exists only if Lambda_tot = 0, so the **boundary
  condition** fixes the total at zero by fiat. This is the Candidate 1 obligation, "quantum asymptotic
  flatness requires a stated physical total-vacuum matching rule" (ledger s08). AdS fixes Lambda < 0
  by its boundary structure. Neither is closure, and neither gives the observed value.

### 2.2 Wheeler-DeWitt: normalizability fixes relations, not Lambda [D; C: C7-C8]

Take k = +1 in conformal gauge, unit-S^3 action (3 pi/(4G)) ∫ d eta (-a'^2 + a^2 - (Lambda/3) a^4),
plus conformal radiation with conserved conformal energy E_r (a state constant). Rescale
x = a sqrt(3 pi/(2G)), lam = 2 G Lambda/(9 pi). The quantum constraint, with one ordering and a
Dirichlet condition at a = 0 (both supplied), is

  -(1/2) psi'' + ((1/2) x^2 - (1/2) lam x^4) psi = E_r psi.

* lam = 0: E_r = 1.5, 3.5, 5.5, 7.5 (exact 2k + 3/2, reproduced to 6e-5).
* lam < 0 (Lambda < 0): normalizable states exist only on discrete curves E_r = E_n(lam). For
  example lam = -0.1 gives 1.653, 4.176, 6.985, 10.02, and lam = -1 gives 2.324, 6.578, 11.65,
  17.32. This is **one relation between Lambda and a state constant**. Lambda would be fixed only if
  E_r were fixed independently, which closure does not do.
* lam > 0 (de Sitter sign, our universe): the potential is unbounded below. The lowest "level" is
  -393 in a box of 12 and -1412 in a box of 16: there is no normalizable spectrum and no
  quantization at all.

### 2.3 Answer to (2)

Gravitational closure (constraint, zero value) makes the origin physical but leaves Lambda_b a
parameter that the constraint ties to the state. Normalizability at most quantizes Lambda relative
to other state constants, and only for Lambda < 0. This does not refute the scoped obstruction. It
confirms Codex's remark that the obstruction is not a no-go for constraint systems, and it shows
that such systems still need an extra rule to fix Lambda.

## 3. QFT in curved spacetime: the same obstruction in local covariant form

### 3.1 Wald's axioms [Q, secondary read; primary not accessed]

The renormalized <T_ab> is required to satisfy: (i) off-diagonal matrix elements between orthogonal
states equal the formal expression; (ii) conservation; (iii) causality / locality; (iv) in Minkowski
spacetime, <0|T_ab|0> = 0. A later fifth condition, no metric derivatives above second order, is
incompatible with the trace anomaly and was dropped. Wald's uniqueness theorem: two prescriptions
differ by a local, conserved curvature tensor. This is Wald 1977 CMP 54, 1; 1978 PRD 17, 1477;
book 1994 pp. 89-91, here as recounted by Tichy-Flanagan §I and Hollands-Wald 2004.

* Massless field: two-parameter freedom, alpha ^(1)H_ab + beta ^(2)H_ab (the variations of
  ∫sqrt(-g) R^2 and ∫sqrt(-g) R_ab R^ab) [Q: HW 2004 essay; Tichy-Flanagan].
* Massive field: Tichy-Flanagan show the Wald-axiom freedom is a free function,
  t_ab = -m^2 R_ab F'(R/m^2) + (1/2) m^4 F(R/m^2) g_ab + ..., with F(0) = 0 required by the
  Minkowski condition. They state that an extra principle would be needed to pin it down [Q,
  §I read].

**Correction to the task framing [D].** Axiom (iv) removes the matter-side c1 g_ab term. It does so
by convention, as a renormalization condition that sets rho_m,vac = 0 in Minkowski. It leaves Lambda_b
in G_ab + Lambda_b g_ab = 8 pi G <T_ab> untouched. It fixes one side of the unphysical split, as
Codex wrote ("predicting one unphysical split does not predict the measurable total").

### 3.2 Hollands-Wald [Q]

Local covariant Wick and time-ordered products, fixed by locality and covariance, scaling,
microlocal spectrum, smoothness and analyticity, are unique up to finitely many parameters at each
order (HW 2001, abstract and intro read). The review's Theorem 2 (HW 2015, pp. 49-51 read) states
that scheme changes are local, covariant, contact-term counterterms D_n of the right dimension.
These include the **c-number terms m^4, m^2 R, R_mu_nu R^mu_nu**. As Lagrangian counterterms they
are a m^4 sqrt(-g), b m^2 sqrt(-g) R, and so on, which shifts T_mu_nu by a' m^4 g_mu_nu + b' m^2 G_mu_nu
+ curvature-squared terms. The 2004 essay (footnote 3, full read) states the additional massive
freedom "m^4 g_ab and m^2 G_ab". The review gives the phi^2 example: the c-number ambiguity
c2 m^2 is fixed "by fixing the expectation value for any given but fixed m = m0, L = L0", which is
an external normalization datum. Stress-tensor conservation for interacting theories is established
in HW 2005 (abstract read). Under the scaling axiom a massless theory has no dimensionful c1, so the
matter-side c1 is a sum of a_i m_i^4 over the massive fields.

### 3.3 Local covariant form of the obstruction [D/A]

Let calA be the HW axioms. If T_mu_nu satisfies calA, so does T_mu_nu + a m^4 g_mu_nu for every
real a. Therefore no rule derived from calA, or from any fixed-background state property of the
resulting theory, fixes a. The observable rho_tot = M^2 Lambda_b + rho_m,vac(a) also contains the
independent gravitational parameter. At order n, D_n adds new constants, so the statement holds at
every perturbative order. This is Codex's obstruction in local covariant form. It is standard, and
the HW 2004 essay states its consequence: "there is still a 'cosmological constant problem'".

### 3.4 Heavy thresholds (Codex Q3) [Q formula (recall); C numbers, retrospective]

The one-loop MS-bar vacuum energy of a real scalar is rho_1 = (m^4/(64 pi^2))(ln(m^2/mu^2) - 3/2).
It runs as d rho/d ln mu^2 = -m^4/(64 pi^2), and the counterterm cancels this running. Compared with
rho_Lambda = 2.52e-47 GeV^4 = (2.24 meV)^4 (Planck 2018 H0 = 67.4, Omega_Lambda = 0.685; computed
rho_Lambda/M^4 = 7.16e-121), m^4/(64 pi^2) is 4.3e30 rho_Lambda for one electron-mass degree of
freedom, 1.5e52 for the Higgs, and 5.6e52 for one top-mass degree of freedom. A rule stated at one
order and one scale must say how these terms are absorbed. Otherwise its value moves by up to 52
orders at the next threshold.

## 4. Unimodular gravity: relocation, not refutation

### 4.1 Trace-free equations [D; C: C5]

Varying only at fixed sqrt(-g) gives M^2 (R_mu_nu - (1/4) R g_mu_nu) = T_mu_nu - (1/4) T g_mu_nu.
The right side is blind to T -> T - C g (C5). Bianchi plus conservation give
d_mu(M^2 R + T) = 0, so M^2 R + T = 4 Lambda_hat, and then M^2 G_mu_nu = T_mu_nu - Lambda_hat g_mu_nu.
Writing T = tau - rho_vac g, the effective constant is

  Lambda_eff = Lambda_hat + rho_vac = (M^2 R + tau^mu_mu)/4,

which is fixed by geometric and non-vacuum data on any slice and is independent of rho_vac.
**The large matter term drops out of the local equations at all orders. The observed constant
becomes an integration constant [D]** (Weinberg 1989 §7: it "does not solve the problem but changes
its nature; why choose this solution" [Q]; Ellis et al. 2011 [Q abstract]).

### 4.2 Henneaux-Teitelboim and Unruh [Q via secondary; C: C4]

S = ∫ [sqrt(-g)(M^2 R/2 - Lambda) + Lambda d_mu calT^mu] + S_m has no background structure. It gives
d_mu Lambda = 0 and d_mu calT^mu = sqrt(-g). T(t) = ∫_Sigma calT^0 d^3x satisfies
T(t2) - T(t1) = (4-volume between the slices), and T is canonically conjugate to Lambda
(Carballo-Rubio et al. 2022, App. A3, eqs. A7-A13). The minisuperspace check gives
d Lambda/dt = 0 and dT/dt = N a^3 (C4). In the quantum theory Lambda plays the role of an energy
conjugate to unimodular time (Unruh 1989; HT 1989; cited, not read), and Delta Lambda Delta V_4 ≳ hbar.

### 4.3 The bridge to §1 [A; structure credited to HT, Unruh, Sorkin]

On a fixed background, c is conjugate to the external time t, which no subsystem reads (DEF 1.1,
1.6). It is a pure phase e^{-ict}. With gravity, rho_vac multiplies the 4-volume, and the phase is
e^{-i rho_vac V_4}. **Whether rho_vac is observable depends on whether the volume is dynamical.**
In GR, sqrt(-g) is dynamical at every point, so rho_vac gravitates and is a free coupling. In
unimodular gravity, sqrt(-g) is fixed (WTDiff) or V_4 is only a global clock (HT). Then rho_vac
decouples locally and survives only as the conserved conjugate of global 4-volume, which is a label
of the state. Unimodular gravity turns "free coupling" into "free conserved charge of the global
state". That is the state-selection slot of FORMAL_CORE §7, where closure supplied nothing (§4 row
"a preferred initial state", BY DEFINITION; COND 5.3; the C1/C2 and CP-12 campaigns).

### 4.4 Radiative stability is contested [Q]; assessment [A]

Padilla (1502.05296 §3.1, full section) argues that the Lagrange multiplier's boundary value absorbs
loop corrections and must be retuned, and Padilla-Saltas (1409.3573, abstract) call classical
unimodular gravity gauge-fixed GR with no new perspective. Carballo-Rubio, Garay and García-Moreno
(2207.08499, §§4-5 and App. A read in part) argue that in WTDiff the fixed volume form protects
Lambda from quantum corrections, making it technically natural. **Assessment:** the dispute is about
which quantity is held fixed. If Lambda_eff is specified by geometric data, no retuning is needed,
but the measured expansion was used as input; Codex: "retuning the counterterm to the already known
expansion is not an explanation". If it is specified as a multiplier boundary value, loops shift it.
In neither case is the value predicted. Codex Q2 fails.

## 5. Weinberg's no-go (Rev. Mod. Phys. 61, 1 (1989))

Assumptions [Q: Padilla §3.2 full section; Weinberg §§3, 6, 7 text read in the authorized Russian
translation, equations not extractable]:

* W1: field content is the metric plus finitely many fields psi_n (any GL(4) tensors).
* W2: field equations follow from a local effective Lagrangian density L(g, psi).
* W3: the target vacuum is translation invariant, with every field constant in spacetime (flat).
* W4: general covariance leaves a residual GL(4) on constant fields, under which L -> det(M) L.
* W5: no fine tuning of parameters.

There are two cases. (a) If the psi equations are independent, then L = sqrt(-g) V(psi), and the
metric equation gives V = 0, which is a tuning. (b) If they are degenerate (Weinberg's assumed
relation 2 g_mu_nu dL/dg_mu_nu = Σ f_i dL/dpsi_i), then L = sqrt(-g) e^{4 phi} V. The equations then
require V = 0 (tuning) or e^phi -> 0, the limit in which all masses vanish. In Weinberg §6, a scalar
sourced by R instead gives G_eff -> 0.

**Escape routes, each negating one assumption:**

* W3 (non-constant fields):
  * a time-dependent scalar on a Minkowski or FLRW vacuum, the Fab Four (1106.2000 abstract: the
    scalar "break[s] Poincaré invariance on the self-tuning vacua");
  * gauge-variant potentials that vary while their field strengths are constant, the 3-form case
    (Weinberg's own §6 footnote [Q translation]);
  * Lorentz-violating vector or tensor backgrounds. Constant tensor vevs alone are still covered by
    the GL(4) argument [A].
* W1 (infinitely many fields): extra dimensions and self-tuning branes, which still meet Weinberg's
  argument below the compactification scale (Padilla §6).
* W2 (nonlocality or global constraints): sequestering, degravitation, Lombriser's global variation,
  and Henneaux-Teitelboim's global time.
* W4 (breaking Diff): WTDiff with a fixed volume form.
* Outside the classical premise: quantum-cosmological measures and selection (Weinberg §§5, 8:
  Hawking 1984, Coleman 1988, the anthropic bound).
* The target: the theorem aims at Minkowski. Self-tuning to an observed nonzero H needs its own
  scale [A].

Weinberg: such theorems rest on "technical" assumptions with interesting exceptions (end of §6 [Q]).

## 6. Conclusions

### 6.1 Is the obstruction correct as scoped?

Yes [D]. On a fixed background every closure property is shift-invariant (§1). With a dynamical
metric the shift is physical, because sqrt(-g) is a dynamical coordinate (§1.5). The local
covariant axioms leave its coefficient free (§3). Gravitational constraints relate it to the state
without fixing it (§2). Unimodular gravity relocates it to the state (§4).

### 6.2 Theorem-like statement: minimal escape class [A, with D components]

Hypotheses:

* **H1.** Matter is a perturbative locally covariant QFT (HW axioms, or Wald's without the Minkowski
  normalization) with finite field content.
* **H2.** The metric couples to renormalized <T_mu_nu> through a local, diffeomorphism-invariant
  effective action Gamma[g] = ∫sqrt(-g)(-rho_tot + M^2 R/2 + O(R^2)) + (loop-generated nonlocal
  terms of light fields). There is no background structure.
* **H3.** All variables are local fields: no global or integrated variables, no measure or state
  selection over parameters, no infinite towers inside the EFT range.
* **H4.** Adjusting fields have constant vacuum values (Weinberg).

Claims:

* **Claim A [D].** Under H1-H2, rho_tot is free. Every value is compatible with the axioms and with
  all fixed-background matter data. It can be fixed only by a renormalization condition (a value
  statement), and the bare parameter must be readjusted at every order and threshold. The shift is
  of order m^4/(64 pi^2) per heavy degree of freedom per unit ln mu^2.
* **Claim B [Q, Weinberg].** Adding H3-H4 adjusting fields cannot reach rho_tot = 0 without tuning,
  except in the limit of vanishing masses.
* **Claim C [A].** A principle P that determines rho_tot without a value statement must contain at
  least one of:
  * **(G)** a global or nonlocal functional sensitive to the constant volume mode: ∫sqrt(-g),
    spacetime averages, global multipliers, or a variable conjugate to global 4-volume;
  * **(B)** background or vacuum-symmetry-breaking structure: a fixed volume form or WTDiff, which
    removes the trace mode from the local equation; or non-constant vacuum fields (a rolling scalar,
    a 3-form potential) that break the vacuum's translation or Poincare invariance;
  * **(S)** Lambda promoted to a dynamical or quantum variable (4-form flux, Henneaux-Teitelboim
    Lambda, wormhole parameters), together with a rule fixing its distribution (boundary condition,
    normalizability, measure, selection);
  * **(X)** structure beyond a finite local 4D EFT: infinite towers or extra dimensions, UV/IR-coupled
    cutoffs, or a fundamental finite Hilbert-space dimension tied to Lambda.

  To give a positive answer to Codex Q2-Q3 it must also contain:
  * **(R)** a scale independent of the matter couplings, such as 4-volume or age, a horizon, a
    history average of the trace, a flux quantum, or ln dim H. The output is
    rho_tot/M^4 = f(scale/M). Every available scale is cosmological, giving the coincidence form
    (H/M)^2 = 3.5e-121 against rho_Lambda/M^4 = 7.2e-121 [C, retrospective]. That form needs H as an
    input, so the target becomes an O(1) number such as Omega_Lambda, which is already exposed.
  * **(P)** protection: the constant mode of <T> enters only through the channel that G or B removes.
    T - (1/4) T g and T - (1/4)<T> g are blind to -C g for every constant (C5-C6). Alternatively the
    value is a conserved or quantized label that local loops do not renormalize (S).

The list is complete relative to H1-H4 **by construction**, since each item negates one hypothesis.
It is not a theorem about quantum gravity beyond those hypotheses. A pure value statement (row R1)
always "escapes" by fiat and has no consequence.

### 6.3 Does GRUT closure supply any ingredient? No [A]

* **G:** "nothing outside" is not a global variational principle. Finite 4-volume is not implied
  (FORMAL_CORE §4, SCOPE row). The closure-fork derivation found that sequestering's finite-volume
  conditions "do not follow from the phrase 'nothing outside'". A global constraint would be a newly
  supplied postulate in the dynamics slot.
* **B:** closure supplies no background structure. HT avoids one only by adding global fields.
* **S:** closure supplies no state selection (FORMAL_CORE §4, §7, COND 5.3). Under unimodular
  gravity or HT, the Lambda question becomes exactly a state-selection question (§4.3). This is the
  most precise GRUT-relevant outcome of this track.
* **X:** the finite d of DEF 1.1 is a regulator. Treating ln d as the de Sitter entropy is the
  Banks-Fischler / Bousso N-bound postulate. It relates Lambda to d and would need
  ln d = 24 pi^2 M^4/rho_Lambda = 3.3e122 [C, retrospective]. Closure does not fix d.
* **R:** closure supplies no scale; constants belong to Theta (FORMAL_CORE §4).
* **P:** closure supplies no protection.
* **The constraint H ≈ 0** comes from diffeomorphism invariance of the supplied dynamics, not from
  closure. In DEF 1.1 time is external and supplied, and the Page-Wootters eps is a convention. Even
  with eps = 0, Lambda remains a parameter (§2).
* **"Every bath internal" in de Sitter:** the horizon bath gives the imported geometric KMS relation
  T = H/(2 pi) (COND 5.2). It relates T to Lambda and does not fix Lambda.

Any GRUT vacuum proposal must therefore add G, B, S or X as an explicit postulate and credit the
matching published class (scorecard), unless it uses an assumption absent from that literature.

## 7. Scorecard (escape class x Codex's four questions)

| Escape class | Q1 extra structure | Q2 fixes residual (sign, scale)? | Q3 stable under heavy thresholds? | Q4 independent consequence / exposed data | Verdict |
|---|---|---|---|---|---|
| R0 fixed-background closure (baseline) | none | no | n/a | none | obstruction holds [D] |
| R1 value statement (Wald axiom iv; asymptotic-flatness BC; Candidate 1 rho_Lambda(mu0) = 0) | renormalization condition or boundary condition | fixes a chosen value (0 conflicts with observation; inserting the observed value is a fit) | only by readjusting the bare term at each order | none | escapes by fiat; not an explanation |
| R2 gravitational constraint closure (compact Sigma, WDW, eps = 0) | reparametrization constraint, zero value | no: relates Lambda to state constants; quantized relations only for Lambda < 0 | n/a | none | origin physical, Lambda free [D/C] |
| R3 unimodular / HT | fixed volume form (WTDiff) or global calT^mu with Lambda conjugate to 4-volume | no: arbitrary integration constant | trace mode decouples at all orders (C5); stability of the residual contested | classically none; quantum Delta Lambda Delta V_4 ≳ hbar | relocates to state selection [D/Q] |
| R4 global / nonlocal (KP sequestering global and local; Lombriser 2019; degravitation) | global variables or averages coupling to ∫sqrt(-g) | global KP: residual set by the history average (1/4)<tau>, small, needs finite 4-volume; local KP: arbitrary constant; Lombriser: Omega_Lambda = 0.704 with a location prior | yes for constant vacuum terms (C6) | KP: transient w ≈ -1 and future collapse (future-testable); curvature claims conditional. Exposed: Planck, DESI BAO, SNe, Omega_K | strongest existing class; credit KP and Lombriser [Q abstracts] |
| R5 Weinberg escapes by non-constant fields or infinite towers (Fab Four, 3-form potentials, self-tuning branes) | rolling or gauge-variant fields, Lorentz-violating vacua, extra dimensions | Fab Four tunes to flat; nonzero residual needs an added scale | unclear; fine tuning tends to reappear; not assessed here | fifth forces, c_T. GW170817 (exposed) strongly constrains Horndeski/Fab Four [recall] | escapes the theorem, not the residual |
| R6 state or measure over Lambda (Baum-Hawking, Coleman, Brown-Teitelboim/Bousso-Polchinski with selection) | Lambda as a quantum variable (4-form flux, wormhole alpha) plus a measure or selection | Hawking/Coleman: Lambda → 0+ (sign and scale fixed, wrong or contested); anthropic: a distribution | flux quanta are not renormalized locally; measure-dependent | anthropic correlations only; no clean held-out test | escapes; a measure is supplied [Q/A] |
| R7 scale-linking (Cohen-Kaplan-Nelson; holographic DE; everpresent Lambda; CosMIn; N-bound / finite dim) | IR scale tied to rho via entropy bounds, discreteness or mode counting | scale ~M^2 H^2 fixed; sign often by construction; O(1) coefficient free except CosMIn's claimed N = 4 pi | usually assumes the large part is removed elsewhere | evolving w(z) or fluctuations; tested only on exposed data (DESI, SNe) | gives the coincidence form; needs H as input [E/Q] |
| R8 GRUT closure as such | none beyond FORMAL_CORE | no | n/a | none | supplies no G, B, S, X, R or P [A] |

## 8. Limits

* **Primary texts not read:**
  * Wald 1977, 1978 and 1994: axioms taken at secondary level from Tichy-Flanagan §I and the HW 2004
    essay. The axiom numbering is not verified.
  * Henneaux-Teitelboim 1989: via Carballo-Rubio et al., App. A3. Elsevier returned 403.
  * Unruh 1989: citation only (APS returned 403).
  * Weinberg 1989: authorized Russian translation, text of §§3, 6, 7 and the start of §8; equations
    were not extractable.
  * Hawking 1984, Baum 1983, Coleman 1988, Brown-Teitelboim, Weinberg 1987, Dvali-Hofmann-Khoury,
    Li 2004: recall only.
  * The one-loop formula: standard, recalled.
  * GW170817 constraints on the Fab Four: recall.
* **WDW minisuperspace:** one operator ordering, a Dirichlet condition at a = 0, finite-difference
  box. Illustrative only. The radiation-energy normalization is defined only in the rescaled variables.
* **Claim C:** complete relative to H1-H4 by construction (argued). Quantum gravity beyond EFT,
  non-perturbative effects and non-Lorentzian UV completions are not covered.
* **Not assessed:** the sign and scale of the sequestering residual, Lombriser's derivation, and
  CosMIn's derivation (abstract level only; a separate `sequestering` track exists).
* **Retrospective inputs:** all numbers using H0 or Omega_Lambda (Planck 2018) are retrospective.
* **Claims:** no new theory and no empirical validation are claimed.

## 9. Sources and read levels

| Source | Identifier | Read level | Used for |
|---|---|---|---|
| Hollands-Wald, local Wick polynomials | CMP 223, 289 (2001), gr-qc/0103074 | abstract; introduction and §3 excerpts on locality/covariance and normal ordering | finitely many parameters; normal ordering not local-covariant |
| Hollands-Wald, stress-tensor conservation | Rev. Math. Phys. 17, 227 (2005), gr-qc/0404074 | abstract | conservation; perturbative agreement |
| Hollands-Wald, "QFT is not merely QM ..." | GRG 36, 2595 (2004), gr-qc/0405082 | full | two-parameter massless freedom; m^4 g, m^2 G for massive; their reframing of the CC problem |
| Hollands-Wald, review | Phys. Rep. 574, 1 (2015), arXiv:1401.2026 | Theorem 2 and examples (pp. 46-51) | c-number counterterms m^4, m^2 R, R^2; fixing by normalization |
| Tichy-Flanagan | PRD 58, 124007 (1998), gr-qc/9807015 | §I | Wald axioms, massive free-function ambiguity, F(0) = 0 |
| Wald | CMP 54, 1 (1977); PRD 17, 1477 (1978); book (1994) §4.6 | not accessed (secondary only) | axioms |
| Weinberg | Rev. Mod. Phys. 61, 1 (1989), DOI 10.1103/RevModPhys.61.1; UFN 158, 639 (1989) (authorized Russian translation) | §§3, 6, 7 and start of §8 text (translation) | no-go, footnotes on exceptions, unimodular, 3-form |
| Padilla, lectures | arXiv:1502.05296 | §§3.1-3.2 full | Weinberg proof; unimodular critique |
| Padilla-Saltas | EPJC (2015), DOI 10.1140/epjc/s10052-015-3767-0, arXiv:1409.3573 | abstract | unimodular gravity is gauge-fixed GR |
| Carballo-Rubio, Garay, García-Moreno | CQG 39, 243001 (2022), arXiv:2207.08499 | §§4-5 (parts), App. A3 | technical naturalness in WTDiff; HT action |
| Ellis, van Elst, Murugan, Uzan | CQG 28, 225007 (2011), arXiv:1008.1196 | abstract | trace-free equations are viable |
| Henneaux-Teitelboim | PLB 222, 195 (1989), DOI 10.1016/0370-2693(89)91251-3 | not accessed (via 2207.08499) | HT action |
| Unruh | PRD 40, 1048 (1989) | not accessed | unimodular time |
| Kaloper-Padilla | PRL 112, 091304 (2014), arXiv:1309.6562 | abstract | sequestering; stated consequences |
| Kaloper-Padilla; KPSZ | arXiv:1406.0711; arXiv:1505.01492 | via the closure-fork packet SOURCES.md | global vs local formulation |
| Lombriser | PLB 797, 134804 (2019), arXiv:1901.08588 | abstract | global variation with respect to M_Pl; Omega_Lambda = 0.704 |
| Charmousis-Copeland-Padilla-Saffin | PRL 108, 051101 (2012), arXiv:1106.2000 | abstract | Fab Four; broken Poincare invariance |
| Bousso-Polchinski | JHEP 0006, 006 (2000), hep-th/0004134 | abstract | four-form discretuum |
| Ahmed-Dodelson-Greene-Sorkin | PRD 69, 103523 (2004), astro-ph/0209274 | abstract | everpresent Lambda |
| Cohen-Kaplan-Nelson | PRL 82, 4971 (1999), hep-th/9803132 | abstract | UV/IR bound |
| Banks | IJMPA 16, 910 (2001), DOI 10.1142/S0217751X01003998, hep-th/0007146 | abstract | cosmological SUSY breaking |
| Bousso, N-bound | JHEP 0011, 038 (2000), hep-th/0010252 | abstract | Lambda linked to Hilbert-space dimension |
| Padmanabhan-Padmanabhan, CosMIn | IJMPD 22, 1342001 (2013), arXiv:1302.3226 | abstract | claimed N = 4 pi determination |
| Hawking; Baum; Coleman | PLB 134, 403 (1984); PLB 133, 185 (1983); NPB 310, 643 (1988) | recall | Lambda → 0 measures |
| Schoen-Yau; Witten | CMP 65, 45 (1979); CMP 80, 381 (1981) | recall | positive ADM energy |

## 10. Files

* `code/obstruction_checks.py`: sections A-D, numpy/scipy/sympy, seed 20261021, 34 checks.
* `results/obstruction_checks.json`: every number quoted above and the check list (34/34 pass).
* Downloads (not in the repository) are under the scratchpad `lit/` and `lit/vac/`.

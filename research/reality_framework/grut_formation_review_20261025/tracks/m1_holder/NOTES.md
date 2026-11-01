# Track m1_holder: Codex's Option 1A curvature phase holder and the isocurvature claims (f8186fa; sextic holder of 94f660d)

1 November 2026 · Claude, issue #13, BRIEF_M1.md · same-session AI work, not outside replication. Packets read only
through `git show f8186fa:…`, `94f660d:…`, `d4747cb:…`, `04c53f5:…`; nothing checked out, nothing committed. **No
newer isocurvature packet exists on the remote**: `git ls-remote origin` and `git log --all --since=2026-10-30` on
2026-11-01 show nothing after f8186fa on any `codex/grut-m1-*` branch, so the brief's "adiabatic residual 4e-12 /
false residual 3.02 / hidden-scalar power 3.435e-70 / 38 checks" could not be read; the analogous quantities were
computed here with our own engine instead (section 4).

**Scope.** One supplied tree-level working model (GRUT-M1 point of d4747cb: H_i = 1.342842735e11 GeV, M_Y = 30 H_i,
h = 40 H_i, μ = 7.159456 H_i (the brief's "7.2" is this number), κ = 3.4587e-3, λ_X = 1e-3, λ_Y = λ_Σ = 1,
ξ_X = −5/6; quadratic inflaton m_I = 3 H_i, Γ_I = 0.5 H_i, prescribed gravity in Codex, Jordan-frame constraint-solved
gravity in ours), homogeneous and classical except for the linear spectator transfer. Nothing here is evidence that
X, Y, Σ exist; nothing is a validation or exclusion. Labels: [D] derived, [C] computed, [E] estimate, [Q] quoted,
[A] argued, [R] recalled.

**Independence.** All numbers below marked "own" come from code written in this track (`code/`): `holder_point.py`
(own potential, own starting guess, BFGS + Levenberg–Marquardt, Hessian by exact sympy differentiation of V(r₁..r₆)
and by a fourth-order finite-difference stencil, 600-start global search, re-derived bound), `rpholder.c` (copy of
our rp_reheating engine `rpreheat.c` with μ → μ_eff(R), initial velocities and full-state restarts; H reconstructed
from the Friedmann constraint, Raychaudhuri H_aux as the geometry check, DP5(4)), `rph.py` (driver, controls,
finite-difference response row, phase slopes, time-shift test), `isocurvature.py` (exact de Sitter Hankel mode
functions with mpmath, own DOP853 check, power table). The only reuse is `rp_reheating/code/rp.py::prepare` (our own
Newton solve of the M1 preparation point). None of Codex's code was run; its numbers are quoted [Q] from REPORT.md,
DERIVATIONS.md, holder_model.json, holder.json, reference.json, perturbations.json, phase_holder.json.

## 1. The operator [D, E]

- V ⊃ (ξR/M_Y)(e^{iα}Σ*Y² + h.c.) = 2 Re[(ξe^{iα}/M_Y) R Σ*Y²]: [R] = 2, [Σ*Y²] = 3, one power of 1/M_Y → the
  operator R Σ*Y² has mass dimension 5 and the term is dimension 4 in the Lagrangian ✓. Writing u(R) = μ + lock·R
  with lock = ξe^{iα}/M_Y and requiring u(12 H_i²) = 10i H_i gives lock = (10i − μ/H_i)/12 in H_i units, i.e. exactly
  Codex's μ_eff(R)/H_i = μ/H_i + (10i − μ/H_i) R/(12 H_i²). Own reconstruction: lock = −0.5966213466 + 0.8333333333 i,
  **ξ = 30.746728738599167, α = 2.192144074179184 rad** — identical to Codex's digits (relative difference 0). u(0) = μ
  (real), u(3 H²) = 5.37 + 2.5i, u(12 H²) = 10i. CONFIRMED.
- The pair (ξ, α) is two real parameters fixed by two conditions (phase −π/2 and |μ_eff| = 10 H_i at R = 12 H_i²
  exactly): Re u(12) = 0 is a **tuned cancellation** of the renormalizable μ by the curvature term at one value of R.
  The equilibrium phase δ_eq = arg u(R) + π rotates from −π/2 (R = 12) through −2.706 (R = 3) to π (R = 0), with
  dδ_eq/d ln R = R Im(lock/u) = 0.7159456 at R = 12 (Codex 0.715945616 ✓). A 1 % change of H² during the last
  e-folds therefore moves the prepared phase by 0.007 rad and (with d ln|N_X|/dδ = 0.554, section 3) the charge by
  0.4 % [E].
- **Size of the coefficient.** At R = 12 H_i² the holder coefficient is |lock|·12 = 12.30 H_i = 1.72 μ: the
  dimension-5 term is larger than the renormalizable cubic of the same field structure. The curvature expansion
  parameter R/M_Y² = 0.0133 is small, so an R² Σ*Y²/M_Y³ term with a similar coefficient would be 1.3 % of the dim-5
  one; but ξR/M_Y² = 0.41 against μ/M_Y = 0.24 shows the curvature series does not have uniform coefficients. One
  loop of Y (vertex μΣ*Y², heat-kernel a₃-type term R·Tr M⁴/M_Y²) generates ξ_loop ≈ c·λ_Y μ/(16π² M_Y) ≈ 1.5e-3·c,
  c = O(1): the required ξ = 31 is ≈ 2e4 times the one-loop size; graviton loops give (M_Y/M_Pl)² = 2.7e-12. So
  ξ ≈ 31 is **not loop-natural**; it must be a tree-level matching input, e.g. a heavy singlet S with
  L ⊃ g₁M_S S R + g₂ S Σ*Y² gives ξ = g₁g₂M_Y/M_S, which needs a large nonminimal coupling (g₁g₂ ≈ 31 for
  M_S ≈ M_Y). NEW [E]; Codex says only "its Wilson coefficient is large" of the sextic and "dimension five … not
  UV-complete" of this one.
- **EFT margin at the prepared point.** |Y|/M_Y = 0.192 (|Y|²/M_Y² = 0.037): the Y amplitude is controlled with a
  factor-5 margin; |Σ|/M_Y = 0.95; but **|X|/M_Y = 2.91** (Codex's 94f660d already noted 2.55 for the old state), and
  the Y-exchange-induced |κ|²|X|⁶/M_Y² is 10 % of λ_X|X|⁴ there, so operators suppressed by M_Y are not
  parametrically small in the X direction. At |Y| ~ M_Y the holder term (≈ 3e5 H⁴) is comparable to M_Y²|Y|² and
  λ_Y|Y|⁴ (8e5 H⁴ each): the fixed-R polynomial stays bounded (section 2) but the truncation is no longer controlled.
  CONFIRMED for Y with the margin stated; caveat for X.
- **CP-odd input.** With κ and μ made real by rephasing Y and Σ, the remaining phase arg(lock/μ) = α is a genuine
  rephasing-invariant CP-odd Lagrangian parameter; it is the supplied sign that our preparation theorem
  (REVIEW_REHEATING.md item 3: no conjugation-even condition can select the sign) requires. Codex states this
  ("the initial phase is encoded in the chosen complex coefficient; this calculation does not derive a unique
  state-selection law"; STABILIZED_CANDIDATE.json scope) and claims no selection law. CONFIRMED.

## 2. Stationary point, Hessian, bound [C]

Own potential V = λ_X|X|⁴ + λ_Y|Y|⁴ + M_Y²|Y|² + λ_Σ(|Σ|² − h²/2)² + 2Re[κY*X³ + u(R)Σ*Y²] + ξ_X R|X|² at R = 12
(ξ_X R = −10), starting guess x² = −ξ_X R/(2λ_X), Y = −κx³/M_Y², Σ = −(u/|u|)h/√2 (not Codex's numbers).

| quantity | own | Codex [Q] | difference |
|---|---|---|---|
| amplitudes (|X|, |Y|, |Σ|)/H_i | 87.19556213287582, 5.754247963285605, 28.38718194436162 | 87.19556213287578, 5.75424796328559, 28.387181944361622 | ≤ 3e-15 |
| arg Σ − 6 arg X; arg Y − 3 arg X | −π/2 (residual 0.0); π | −π/2 | exact |
| gradient norm | 7.3e-12 (exact), 1.1e-7 (FD of V) | 1.3e-12 | — |
| V at the point / H_i⁴ | −32480.58531140619 | −32480.585311406205 | 4e-16 |
| Hessian eigenvalues (sympy) | 7e-14, 13.091408829452917, 14.589249978667548, 538.1466434275553, 1546.659721973108, 3239.8852511586965 | 4.7e-14, 13.091408829453263, 14.589249978666821, 538.1466434275541, 1546.6597219731073, 3239.885251158698 | ≤ 5e-14 |
| Hessian eigenvalues (4th-order FD) | −1.3e-6, 13.0914109, 14.5892531, 538.146641, 1546.659714, 3239.885249 | | ≤ 2e-7 |
| gauge null vector T·r | H·(T r)/|T r| = 1.1e-13 | "remaining eigenvalue is the gauge direction" | ✓ |
| 1 − 2η² dV/dR (effective Planck ratio) | 1 + 4.80667e-11 (dV/dR = −7902.45) | 1 + 4.80667e-11 | ✓ |

Hand check at δ = −π/2 (X real, Y real negative, Σ = −is): 2λ_X x³ − 3κyx² + ξ_X·12·x = 0, (M² + 2λ_Y y²)y = κx³ +
2|u|sy, 2λ_Σ(s² − h²/2)s = |u|y², residuals ≤ 2e-11; the point is a local minimum (all five physical eigenvalues
positive). Mode content: 13.09 = X radial (0.99 ReX, 0.15 ReY), 14.59 = X phase mixed with Σ (0.89 ImX, −0.46 ReΣ;
the gauge mode is the orthogonal combination), 538 = Y radial, 1547 = Y phase, 3240 = Σ radial.

**Bound.** Codex's chain re-derived: reserve A = 0.2λ_X; min_x[0.8λ_X x⁴ − 2|κ|yx³] = −27|κ|⁴y⁴/(16(0.8λ_X)³) (own
minimisation; the inequality was also tested on 2e4 random (x, y) with no violation), so B = λ_Y − 0.47165 =
0.5283522987 > 0; then By⁴ − 2|u|sy² ≥ −|u|²s²/B, Ax⁴ − a_x x² ≥ −a_x²/(4A), and the Σ quartic absorbs the rest:
V ≥ λ_Σh⁴/4 − a_x²/(4A) − (λ_Σh² + |u|²/B)²/(4λ_Σ) = 640000 − 125000 − 800369.7 = **−285369.684 H_i⁴** (Codex
−285369.684 ✓). It holds for every fixed R (own values: −79965 at R = 0, −62035 at R = 3, −96554 at R = −1.063, the
lowest R reached in the run). It is loose: a 600-start L-BFGS-B search in the six real fields (log-uniform radii to
200 H_i, random phases) finds **one** distinct minimum, the prepared point itself (V = −32480.585), and 2e5 random
points never go below −25987. So the prepared point is the global minimum at R = 12 to the resolution of the search,
8.8× above the bound. No counter-direction: the quartic form is positive after the mixed-quartic inequality and
dominates the cubic (|u|s y²) and quadratic (a_x x²) terms. The F > 0 caveat for huge fields (F = M_Pl² − 2W needs
s y² ≳ 8e13 H_i³) is outside any EFT regime, as Codex says. CONFIRMED.

**Sextic holder of 94f660d** (ξ₆R(iΣ*X⁶ + h.c.)/Λ⁵, Λ = 1e15 GeV, ξ₆ = 4e8, c₆ = 2.0959e-10 H_i⁻⁵): own stationary
point (86.6572, 3.8228, 28.3376) H_i, δ = −2.099927, lightest physical eigenvalue 10.66925 H_i², holder term
−4342 H_i⁴ at the point; Codex (86.657, 3.823, 28.338), −2.09993, 10.669 ✓. CONFIRMED (and it does not reach −π/2,
as Codex says; Option 1A does, by the tuned cancellation of section 1).

## 3. Charge change with the holder, own engine [C]

`rpholder.c` = our `rpreheat.c` with u = μ + lock·R, where R is the engine's own trace-equation Ricci scalar
(signature −+++, R = 12H² in de Sitter, R = 2m_I²I² − J² + O(η²) in the oscillation era, 0 in radiation; this is the
same quantity and sign as Codex's prescribed 18I² − J², checked: own fine-sampled minimum R = −1.0636 at τ = 0.70
against Codex's range −1.0633 … 12). The holder's R-dependence is treated as an explicit time dependence of V (its
contribution to F is dropped, relative size 2η²|lock s y²| ≈ 1e-11; Codex drops all backreaction, η² = 3.0e-15).
Initial states: original = our own Newton solve of the M1 preparation (76.39978266, −2.11150985 − 0.93932904i,
−28.29314213i; identical to d4747cb), holder = own stationary point of section 2; zero velocities, a = 1, H = 1,
I₀ = √6/3 M_r from the constraint; rtol 1e-10/atol 1e-13; late mean τ-weighted over the last 20 % (W80) as Codex.

| run | own N_X (W80, 640–800) | Codex [Q] | rel. diff | own W90 (720–800) | Codex W90 |
|---|---|---|---|---|---|
| original (lock = 0) | −8120.740082 | −8120.739696 | 4.8e-8 | −8118.614873 | −8118.614541 |
| holder | −8267.437810 | −8267.437032 | 9.4e-8 | −8267.128260 | −8267.128654 |
| **ratio** | **1.01806458 (+1.806458 %)** | **1.01806453 (+1.806453 %)** | 4.7e-8 | 1.018293 | 1.018295 |

Controls: rtol 1e-12 leaves the ratio at 1.018064576 (change 4e-13); Gauss ≤ 4.7e-12, |H_aux − H|/H ≤ 5.4e-9, hidden
work ledger (now including dV/dR·Ṙ with the prescribed Ṙ = 6m_I²IJ + 2(3H + Γ)J²) ≤ 3.5e-11 (original 2.7e-14).
H_i refit at fixed mass per charge, H_new/H_old = ratio^{−2/3}: own **0.98813538 → 1.32691042e11 GeV** (Codex
0.98813541 → 1.32691046e11 ✓); fully-cooled hypothesis 0.98263381 (Codex 0.98263384 ✓). CONFIRMED.

Two things Codex's packet does not say:
- **Window dependence is 1e-3, not 3.7e-5.** Codex's 3.7e-5 is the W80→W90 change of the holder N_X alone; the
  *ratio* moves from 1.01806 (640–800) to 1.01829 (720–800) and to **1.01961 (800–1000)** (own end-1000 runs:
  holder −8277.72 vs original −8118.50). The holder run's |N_X| is still growing at τ ~ 1000 (+0.12 % per 200 H_i t)
  while the original's falls (−0.03 %). The source is physical in this model: the holder run leaves **N_YL = 3.08
  light-Y free quanta** (original 0.0068; same diagnostic as rp_reheating) whose charge-3 content (≈ 9 units ≈ 0.11 %
  of N_X) is slowly exchanged with X through κ. So "+1.806 %" is the 640–800 window value with an intrinsic
  ±0.15 % late drift, and the H_i refit carries ±0.1 % from this alone. CORRECTED (precision statement).
- **The holder excites Y by 450× [C, E].** The π/2 rotation of arg μ_eff as R falls from 12 to 0 within ~1/H_i is a
  non-adiabatic kick to the Y condensate (ω_YL = 22 H_i): N_YL 0.0068 → 3.08, ε −0.199 → −0.152. Under the
  rp_reheating/hw_verify debris ledger scaling (m_I = 10 rows: N_YL 1.6–1.7 ↔ R_central 5–8) this is R_central ≈ 10
  [E] against the energy-conserving proxy's 7–8 and the instant release's 1.14: the Option 1A holder moves the quiet
  (3, 0.5) point (R ≈ 0.008 in rp_reheating) into the marginal class under that proxy. Not computed with the full
  ledger here; flagged as NEW for the debris track. Charge build-up: holder N_X(τ)/N_X(800) = 0.62, 0.71, 0.92, 1.01,
  1.07, 1.06, 1.00 at τ = 0.5, 1, 2, 3, 5, 10, 20 (original 0.78, 0.87, 0.91, 1.02, 1.02, 1.02, 1.00).

**Response row and phase slope.** Own central differences (steps 1e-2, 1e-3, Richardson) of ln|N_X| along the five
unit eigenvectors of our Hessian, position at zero velocity and velocity at zero displacement, in canonical real
fields r = √2(Re, Im): (−0.0154820, −0.0301824, 0.00187089, −0.000376023, −7.82e-6; −0.00635517, 0.00707972,
−6.1500e-5, 6.804e-7, −3.86e-8). Codex's canonical 12-row [Q perturbations.json] projected on our eigenvectors gives
(−0.0154820, −0.0301825, 0.00187090, −0.000376028, −7.82e-6; −0.00635519, 0.00707975, −6.1501e-5, 6.806e-7,
−3.89e-8): norm relative difference **3.9e-6** (the 1e-3 and 1e-2 steps agree to the same level; the two smallest
entries differ by 2e-4 and 5e-3 relative, within finite-difference noise on 1e-7-sized numbers). Gauge direction:
1.8e-10 (Codex 7e-12). Σ-phase rotation at fixed amplitudes and velocities: own d ln|N_X|/dδ = 0.5543820 (steps
1e-3, 5e-4) vs Codex 0.5543843 (tangent and finite differences): 4e-6. For the original history the same
Σ-only rotation gives −0.11650, while re-preparing the conditional minimum at the shifted held phase (Y re-solved,
the e160356/d4747cb definition) gives −0.2199687 and −0.2199687 (steps 1e-3, 5e-4) = d4747cb's −0.21996870 /
−0.21996867 exactly; the two definitions are different perturbations and both reproduce. CONFIRMED.

## 4. Isocurvature [D, C, A]

### 4a. A pure time shift is the adiabatic mode and carries no charge-to-entropy perturbation at k → 0 [D]

Let the complete background system (inflaton I, J, radiation ρ_r, scale factor, the six real scalars and their
velocities, with R = R(I, J, ρ_r, …) and H = H(…)) be written as an autonomous first-order system ẏ = f(y), as it is
in Codex's evolve.c and in our engine: nothing in f depends on t explicitly, all time dependence of μ_eff enters
through R(y). If y(t) solves it, so does y(t + ε) for every constant ε; expanding, δy(t) = ε ẏ(t) = ε f(y(t)) satisfies
the linearised (tangent) equation δẏ = Df(y) δy identically, *including* the term −V_rR δR with δR = ε Ṙ — the
"curvature forcing". Any scalar function g(y) of the state then has δg = ε ġ, so its fractional perturbation is
δg/g = ε d ln g/dt. For two comoving densities n_A, n_B (the X charge N_X = −2a³Im(X*Ẋ) and the comoving radiation
entropy S ∝ (ρ_r a⁴)^{3/4}) the entropy perturbation S_AB ≡ δn_A/n_A − δn_B/n_B = ε (d ln n_A/dt − d ln n_B/dt)
vanishes once both are separately conserved, which is the late state of every reheating run (N_X constant up to the
slow κ exchange, comoving entropy constant after Γ_I decay); on any slicing in which the two are compared at the
same local time it vanishes identically. In general relativity at k → 0 this is the first half of Weinberg's theorem
(astro-ph/0302326): the time-shift (plus spatial rescaling) of the background is always a solution of the
perturbation equations, with ζ and R constant and all species sharing δρ_A/ρ̇_A = δρ_B/ρ̇_B; the gauge-invariant
form S_AB = 3(ζ_A − ζ_B) with ζ_A = −ψ − Hδρ_A/ρ̇_A is Wands–Malik–Lyth–Liddle (astro-ph/0003278). The second half
of the theorem (that this mode is the k → 0 limit of a physical finite-k mode) needs the metric constraints, which
the 12×12 engine does not contain — Codex says so.

Consequence for Codex's engine: feeding δχ = εχ̇ into the tangent system *without* δR = εṘ is not a solution, so a
non-zero "residual" there measures the forcing term, not an isocurvature mode; with the forcing it must vanish to
integration accuracy. The brief's quoted 4e-12 (with forcing) and 3.02 (without) have this expected structure; the
packet itself is not on the remote, so the numbers are NOT CHECKED. Own nonlinear version (`rph.py`): take the holder
solution at t₀ and t₀ + ε and restart to τ = 800.

| t₀, ε | full shift d ln N_X/dε | charge − entropy, full | scalar-only shift | background-only shift | sum |
|---|---|---|---|---|---|
| 3, 1e-2 | −8.95e-7 | −8.95e-7 | −0.4789 | +0.4812 | 2.3e-3 |
| 3, 1e-3 | −8.95e-7 | −8.95e-7 | −0.4806 | +0.4809 | 2.3e-4 |
| 6, 1e-2 / 1e-3 | −8.99e-7 | −8.99e-7 | −0.2475 / −0.2480 | +0.2481 / +0.2480 | 5.6e-4 / 5.6e-5 |
| 12, 1e-2 / 1e-3 | −9.03e-7 | −9.03e-7 | −0.1269 / −0.1271 | +0.1271 / +0.1271 | 1.6e-4 / 1.5e-5 |

The full shift gives −9e-7 per unit H_i t at every t₀ and ε — exactly the late drift rate d ln N_X/dt of the window
(the W80→W90 change 3.7e-5 over 80 H_i t = 4.7e-7), i.e. the expected numerical zero; δ ln S is 1e-10. The
scalar-only shift is a charge–entropy perturbation of order H(t₀)·ε (0.48, 0.25, 0.13 ∝ 1/t₀), the
background-only shift is its negative, and their sum vanishes linearly in ε. CONFIRMED (as a theorem and
numerically).

### 4b. de Sitter suppression of the heavy-spectator power [C]

Own exact mode functions: for constant m in de Sitter, S(x) = √(π/2) e^{−πμ/2} x^{3/2} H⁽¹⁾_{iμ}(x), μ = √(m² − 9/4),
x = k/(aH) = e^{−N}, S_N = −x dS/dx (mpmath, 40 digits; massless check |S|² → 1.0000000021; own DOP853 integration
from x₀ = 1000 and 2000 agrees to 1.8e-8). For m² ≥ 13 H² all five modes are over-damped oscillators: |S|² ∝ e^{−3N}
with a residual oscillation of relative amplitude ~e^{−πμ} = 3.2e-5 and 1.6e-5 for the two light modes (μ = 3.29,
3.51) and < 1e-31 for the rest. Power P(N) = Σ_a |L_pos,a S_a + L_vel,a S_N,a|²/(4π²):

| N | own L (own modes) | Codex's L on our basis (own modes) | Codex table [Q] | ratio/e^{3ΔN} from previous row (Codex / own) |
|---|---|---|---|---|
| 0 | 1.985802e-5 | 1.985817e-5 | 1.985817e-5 | — |
| 2 | 4.984098e-8 | 4.984135e-8 | 4.984135e-8 | 0.98760 / 0.98760 |
| 4 | 1.235743e-10 | 1.235752e-10 | 1.235752e-10 | 0.99975 / 0.99975 |
| 8 | 7.592291e-16 | 7.592348e-16 | 7.592348e-16 | 1.00005 / 1.00005 |
| 20 | 1.761177e-31 | 1.761190e-31 | 1.761190e-31 | 0.99993 / 0.99993 |
| 50 | 1.443003e-70 | 1.443013e-70 | 1.443013e-70 | 1.00007 / 1.00007 |
| 60 | 1.350366e-83 | 1.350376e-83 | 1.350376e-83 | 0.99996 / 0.99996 |

With Codex's response row the table reproduces to 1.4e-9 at every N (their ODE + analytic continuation agrees with the
exact Hankel functions); with our finite-difference row to 7.5e-6. The e^{−3ΔN} law holds to 1e-4 for N ≥ 2 (the
0→2 step still feels the gradient term); the 1e-4 deviations are the e^{−πμ} oscillations above. √P(50) = 1.2e-35.
CONFIRMED. The "3.435e-70" of the unpushed packet is 2.38× the f8186fa value and matches neither an incoherent
position+velocity sum (1.27e-70) nor a position-only convention (6.0e-71) of the same row: NOT CHECKED.

### 4c. What remains for a real CMB isocurvature bound, and is "no resolved signal" expected? [A]

Expected, yes, on general grounds: every scalar direction at the prepared point has m² ≥ 13 H_i² > 9H²/4 during
inflation, so there is no light spectator; superhorizon fluctuations of all six fields decay as a^{−3/2} and the
field configuration becomes a deterministic function of the local background (R(t)) — the heavy sector tracks the
local clock. By 4a that is the adiabatic mode, and the only intrinsic charge–entropy perturbation is the exponentially
suppressed remnant of 4b (e^{−150} for 50 e-folds). Codex's conclusion is therefore the generic expectation for a
heavy spectator (Weinberg's theorem plus a^{−3} decay), not a feature of this model. What a real P_CDI needs beyond
this: (i) the response of the *final* dark-matter density per entropy to a local charge perturbation — the mass per
charge of the formed objects, the fractions η, p and the gas/core split as functions of local N_X/s (the formation
transfer function Codex names "response of the final mass per charge"); (ii) the metric, inflaton and radiation
perturbations and δR in the transfer, so that the observable is the gauge-invariant S_{Xγ} = 3(ζ_X − ζ_γ) rather than
the fixed-slicing δN_X/N_X (the adiabatic subtraction of 4a, which the 12×12 engine cannot perform); (iii) the vector
sector: with e = 1e-6 the gauge-boson mass at the prepared point is m_A = e√(2Σq²|z|²) ≈ 2.7e-4 H_i, so the
longitudinal/Goldstone direction is effectively massless during inflation and acquires an H/2π spectrum per e-fold —
pure gauge at k = 0 (Codex's Ward check 7e-12) but a physical longitudinal mode at finite k whose late charge-density
imprint is not computed; (iv) the inflationary history itself (50 e-folds at constant H are assumed; δ_eq drifts
0.72 rad per e-fold of ln R, section 1, so slow-roll variation of H² during the last e-folds modulates the prepared
phase — adiabatic, but a sizeable correction to the homogeneous charge); (v) the correlated (adiabatic) part — the
charge yield's dependence on the local clock is adiabatic, but a scale-dependent transfer (the +9.9 %/+1.8 % of the
two holders show the yield is history-sensitive) must still be shown to reduce to ζ at k → 0 through the formation
stage; (vi) small-scale modes leaving the horizon in the last 2–4 e-folds keep fractional charge power 1e-8–1e-10
(table above), irrelevant for the CMB but a seed for inhomogeneous formation at comoving scales ~e^{−46} of CMB
scales. Verdict: the conclusion "no resolved correlated charge-entropy signal" is CONFIRMED as an expectation
[A], not as a computed bound; Codex states the same limits.

## 5. Verdicts

| # | Claim (Codex) | Label | Codex [Q] | Own | Verdict |
|---|---|---|---|---|---|
| 1 | Operator dimension 5; μ_eff(R) ⇔ (ξ, α) | D | ξ = 30.74672874, α = 2.19214407 | identical (0 difference); μ = 7.159456 H_i | CONFIRMED |
| 2 | Size of ξ ≈ 31 | E | "dimension five, not UV-complete" | 2e4 × one-loop Y-loop size; tree-level input (g₁g₂M_Y/M_S); tuned Re u(12) = 0; dim-5 term = 1.72 μ | NEW |
| 3 | EFT control at the prepared point | E | not stated | |Y|/M_Y = 0.19 (margin 5), |Σ|/M_Y = 0.95, |X|/M_Y = 2.91; κ²X⁶/M_Y² = 0.10 λ_X X⁴ | CONFIRMED for Y, caveat X |
| 4 | Phase encoded in the complex coefficient; no selection law | D | stated | α rephasing-invariant CP phase; consistent with REVIEW_REHEATING item 3 | CONFIRMED |
| 5 | Stationary point (87.19556, 5.75425, 28.38718), δ = −π/2 | C | as stated | same to 3e-15; δ = −π/2 exactly; hand equations ≤ 2e-11 | CONFIRMED |
| 6 | Hessian (13.0914, 14.5892, 538.147, 1546.66, 3239.89) + gauge zero | C | as stated | sympy exact: same to 5e-14; FD 2e-7; null vector T r | CONFIRMED |
| 7 | V ≥ −285370 H_i⁴ at R = 12; bounded for all fixed R | C/D | −285369.684 | chain re-derived, −285369.684; global search: prepared point is the only minimum (−32480.585) | CONFIRMED (loose by 8.8×) |
| 8 | Effective Planck ratio 4.81e-11 | C | 1 + 4.80667e-11 | 1 + 4.80667e-11 | CONFIRMED |
| 9 | Sextic holder 94f660d: δ = −2.09993, lightest 10.669 | C | as stated | −2.099927, 10.66925, amplitudes to 4 digits | CONFIRMED |
| 10 | Charge −8267.437 vs −8120.740, +1.80645 % | C | 1.01806453 | own engine 1.01806458 (W80), rtol 1e-12 unchanged | CONFIRMED |
| 11 | Window sensitivity "3.7e-5" | C | 3.7e-5 (N_X alone) | ratio: 1.01829 (W90), 1.01961 (800–1000); N_YL = 3.08 vs 0.0068 explains the drift | CORRECTED (ratio ±1.5e-3) |
| 12 | H_i refit 0.98813541 → 1.32691046e11 GeV | C | as stated | 0.98813538 → 1.32691042e11 (±0.1 % from 11) | CONFIRMED |
| 13 | Phase tangent 0.554384 | C | 0.5543843 | 0.5543820 (Σ-only); original −0.1165 Σ-only / −0.2199687 re-prepared (= d4747cb) | CONFIRMED |
| 14 | Response row (perturbations.json) | C | 12-row | own FD row on the same eigenvectors: 3.9e-6 norm difference; gauge 1.8e-10 | CONFIRMED |
| 15 | Finite-k controls 7.33e-5, 7.33e-7 | C | as stated | not run | NOT CHECKED |
| 16 | Time shift = adiabatic mode, zero charge–entropy perturbation | D/C | 4e-12 / 3.02 (unpushed) | theorem (section 4a); own: full shift −9e-7 (= ε Ṅ/N), scalar-only −0.48/−0.25/−0.13 at t₀ = 3/6/12 | CONFIRMED (numbers of the new packet NOT CHECKED) |
| 17 | Power 1.236e-10 (4), 7.59e-16 (8), 1.76e-31 (20), 1.44e-70 (50) | C | as stated | exact Hankel: 1.4e-9 with Codex's row, 7.5e-6 with own row; e^{−3ΔN} to 1e-4 | CONFIRMED |
| 18 | Hidden-scalar power 3.435e-70 (unpushed) | Q | 3.435e-70 | not on the remote; 2.38× f8186fa; matches no tested convention | NOT CHECKED |
| 19 | "No resolved correlated charge–entropy signal"; proxies, not P_CDI | A | as stated | expected for a heavy spectator; missing items (i)–(vi) listed | CONFIRMED as expectation |
| 20 | Holder's Y excitation | C/E | not reported | N_YL 0.0068 → 3.08; R_central ≈ 10 under the rp ledger scaling [E] | NEW |
| 21 | Gas, BGK hierarchy, CLASS, FOLLOWUP record-selection | — | — | outside this track | NOT CHECKED |

## 6. Limits

- Homogeneous classical background with a quadratic inflaton surrogate; the holder's contribution to F and all
  backreaction are O(η²) = 3e-15 and dropped or kept only through the Friedmann constraint; no finite-k gauge
  modes, no metric perturbations, no decays, no thermal corrections — the same scope as Codex's calculation.
- The response row and the power table use the fixed-slicing δN_X/N_X and Codex's normalisation
  P = Σ|L·S|²/(4π²); they are spectator proxies, not P_CDI (section 4c).
- The global-minimum statement is numerical (600 starts, radii ≤ 200 H_i), not a proof; the proof is Codex's bound.
- The loop and tree estimates of section 1 are order-of-magnitude [E]; no matching calculation.
- The debris estimate R_central ≈ 10 for the holder's N_YL = 3.08 reuses the rp_reheating/hw_verify ledger scaling
  without re-running the ledger [E].
- Same-session AI work with own code; not outside replication.

## Reproduction (≈ 2 min on one core)

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
T=research/reality_framework/grut_formation_review_20261025/tracks/m1_holder/code; cd $T
python3 -B holder_point.py        # parts 1-2 -> ../results/holder_point.json, stationary_point.json   (9 s)
python3 -B rph.py                 # part 3 and 4a, ~110 C runs -> ../results/charge_transfer.json     (80 s)
python3 -B isocurvature.py        # part 4b -> ../results/isocurvature.json                            (5 s)
```
Credit: Weinberg 2003 (astro-ph/0302326) and Wands–Malik–Lyth–Liddle 2000 (astro-ph/0003278) for the adiabatic-mode
theorem; Bunch–Davies de Sitter modes and Hankel-function solutions (standard, e.g. Birrell–Davies); Dine–Randall–Thomas
Hubble-induced terms [R]; Codex's bound chain (DERIVATIONS.md) re-derived, not reused.

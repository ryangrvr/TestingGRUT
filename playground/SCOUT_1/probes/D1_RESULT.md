> **W4-0 (owner ruling): TC-4 status → HOSTILE TEST OPEN — OPERATOR-GENERATION LOOPHOLE.** Not promoted while
> W4-OG is open. Wilsonian RG generically generates every symmetry-allowed operator. A zero bare coefficient of an
> allowed operator may therefore be a *tuning*, not a supplied operator-class choice. The five instances are kept as
> evidence only.

# SCOUT-1 D1-SCOUT RESULT — does inherited GRUT structure select a nonlinear universality class?

**SCOUT-ONLY premise experiment** (owner ruling after ZOOM_OUT_04). Not a GRUT premise adoption. GRUT-RAI is
untouched.

**Charter:** `PROBE_CHARTERS.md` §D1-SCOUT.

**Files:**
- `d1_scout.py`;
- logs: `d1_scout_A.log` (exact transfer operator), `d1_scout_BC.log` (operator scan + pin role), `d1_scout_cubic.log`
  (stable rerun of the cubic current).

**Labels:**
- **B. OPERATOR-CLASS-SUPPLIED** at the level of the D1 hypothesis class;
- **C. NONLINEARITY-FORBIDDEN** for the KPZ operator within inherited symmetry;
- KNOWN RESULT IMPORT: KPZ / driven diffusive systems; Model A/B; 1D φ⁴ transfer operator — STANDARD ✓;
- theorem candidate **TC-4 (selector relocation no-go)**.

## 0. Verdict

> **Not A.** Inherited GRUT structure does not force a nonlinear universality class.
>
> **(C) Within inherited structure, the relevant KPZ operator is forbidden or absent.**
> - The S2 drift `−Kx − 4βx³` is gradient, odd (Z₂: x → −x) and on-site. A current term `λ₂u²` is Z₂-even, so it is
>   forbidden in the equation for a Z₂-odd field.
> - Reciprocal K (E-1) gives no drive at all.
> - The affinity ring (E-4) supplies a drive, but with constant mobility the current is linear (a Galilean shift).
> - The inherited operator content therefore stays in the Gaussian/EW class or is gapped.
>
> **(B) Within the D1 hypothesis class, the IR class is decided by which operator is supplied.** `λ₂u²` gives KPZ;
> `λ₃u³` gives the Z₂-odd cubic class (see §1); no nonlinear current gives EW. Several equally admissible operators
> lead to different IR classes.
>
> **The pin reveals; it does not select.** Under identical λ₂ dynamics, pin > 0 cuts the scaling off: W saturates at
> 2.7 (r = 0.01) and 1.3 (r = 0.1), against 9.8 and still growing at r = 0. With the inherited S2 quartic, pin = 0
> produces **no critical IR at all**: the quartic regenerates a mass (§1A).

## 1. Results

**A. Inherited quartic on a pin-free chain** (exact 1D transfer operator, T = 1):

| β | 1 | 0.1 | 0.01 | 10⁻³ | 10⁻⁴ |
|---|---|---|---|---|---|
| ξ | 0.701 | 1.326 | 2.731 | 5.815 | 12.494 |
| ξ·β^{1/3} | 0.701 | 0.615 | 0.588 | 0.582 | **0.580** |

So `ξ ≈ 0.58 β^{−1/3}`: finite for every β > 0. Only β = 0 is gapless (an EW interface with a zero mode).

**B. Conserved scalar with bond noise, operator scan** (N = 4096, T = 1; growth exponent of the integrated height
over t windows [5–50, 50–500, 500–2000]):

| Operator content | Origin | Exponents | Class |
|---|---|---|---|
| J = 0 | E-1 reciprocal (inherited) | 0.243, 0.241, 0.277 | EW (1/4) |
| J = εu | E-4 affinity, constant mobility (inherited) | 0.241, 0.260, 0.248 | EW |
| Model B, S2 quartic made conservative | inherited operator, conserved form (D1) | 0.250, 0.245, 0.214 | EW / diffusive |
| affinity + Model-B quartic | inherited + D1 | 0.267, 0.233, 0.233 | EW |
| J = εu² | **supplied** density-dependent mobility (breaks Z₂) | 0.326, 0.298, 0.392 | KPZ (1/3; noisy) |
| J = εu³ | **supplied** Z₂-odd cubic current | first run unstable (centred flux). Lax–Friedrichs rerun (ε = 0.2): 0.242, 0.241, 0.239 | EW-like (marginal operator; log corrections not resolved) |

**B′. Cubic rerun** (`d1_scout_cubic.log`, Lax–Friedrichs flux, ε = 0.2):
- The cubic current gives 0.242 / 0.241 / 0.239, i.e. the EW-like marginal class.
- The λ₂ reference with the same flux and ε = 0.2 gives 0.256 / 0.250 / 0.277. It is still crossing over by
  t = 2000, as expected from `t_x ∝ ε⁻⁴` (W3-NL) plus the extra numerical diffusion of the Lax–Friedrichs scheme.
- So even inside the D1 class, **which** IR class is visible on a finite window depends on supplied magnitudes
  through the crossover scale.

**C. Pin role (λ₂ dynamics):**

| pin r | W(t = 10, 100, 1000, 2000) |
|---|---|
| 0 | 1.71, 3.68, 7.39, 9.79 |
| 0.01 | 1.67, 2.59, 2.74, 2.73 |
| 0.1 | 1.28, 1.32, 1.30, 1.34 |

## 2. D1-3 map: assumption → operator → IR class

| Assumption | Status | Operator it allows / forbids | IR class (1D) |
|---|---|---|---|
| reciprocity / no drive (E-1) | earned (E-1 passivity) | forbids all drive currents | EW / diffusive |
| cycle affinity (E-4) | declared class (supplied model) | allows a drive. **Linear** without a mobility nonlinearity | EW (Galilean shift) |
| Z₂ x → −x (S2 odd drift) | inherited from the declared drift form (supplied choice) | forbids `λ₂u²`; allows `λ₃u³` | EW, or cubic (marginal) |
| density-dependent mobility / exclusion | **D1-only (newly supplied)** | `λ₂u²` when Z₂ is broken | KPZ |
| conserved vs site noise | site noise in the record (S-8); conserved noise is D1-only | conserved density vs non-conserved order parameter | Model B vs Model A |
| on-site quartic (S2) | inherited | mass regeneration (non-conserved); irrelevant (conserved, diffusive) | gapped / EW |
| detailed balance | uniform `T_i` only (supplied profile; FDT borrowed) | gradient dynamics forbids non-potential terms | EW |
| range of couplings / noise | supplied (graph S-1) | long-range ⇒ continuous exponents (W3-TC3) | continuum |
| dimension | supplied | φ⁴ critical point exists only for d ≥ 2 (tuned) | Ising Model A/B at tuned r_c |
| momentum conservation | absent (overdamped record) | — | — |
| boundary class | supplied (W1-S) | edge exponents | — |

## 3. D1-5 information accounting

| Item | Before flow | After flow |
|---|---|---|
| UV numerical couplings (ε, β, D, T magnitudes) | free | **forgotten** in the IR exponent. They survive as crossover scales and amplitudes (G-moved) |
| operator content (λ₂ vs λ₃ vs none) | — | **still supplied**; it decides the class |
| symmetry class (Z₂, parity, detailed balance) | — | **still supplied** |
| noise class (conserved / site) | — | **still supplied** |
| dimension | — | **still supplied** |
| universal IR number | — | fixed **given** the above (1/4 or 1/3) |

**Not a TRUE NET SELECTION.** Specifying the IR law still needs the operator class. The flow erased magnitudes, not
choices. **D1 exchanged a supplied IR label (z, β) for a supplied operator / symmetry class.**

## 4. D1-6 GRUT bridge

The only admitted nonlinearity, S2's `−4βx³`, extends locally to φ⁴ dynamics: Model A for site noise, Model B for
conserved noise. Its consequences:
- **1D:** gapped (non-conserved) or diffusive EW (conserved).
- **d ≥ 2:** the Ising class, but only at a **tuned** `r_c(β, T) ≠ 0`. The tuning moves from pin = 0 to `r_c`.

The KPZ operator **cannot** be obtained by extending S2 under locality / conservation and the inherited Z₂. It must be
**newly supplied** (Z₂ breaking + drive + a density-dependent mobility). **This is not GRUT-internal selection.**

## 5. Theorem candidate TC-4 — SELECTOR RELOCATION NO-GO (candidate)

> **RG universality can erase the magnitude of a relevant coupling (it survives only as a crossover scale, a
> G-moved datum, TC-1), but it cannot select which relevant operator / symmetry class is present. The IR class is a
> function of upstream operator content.**

**Instances:**
- W3-NL: ε forgotten; the asymmetric operator is supplied.
- D1: inherited Z₂ forbids λ₂; adding it is a supply.
- W1-R: the member m is fixed by the symmetry class.
- W2-ETH: the RMT class is fixed by the antiunitary class.
- W2-DT: g is traded for a scale; the field content (asymptotic freedom) is supplied.

**Grade:**
- CANDIDATE: a classification over five instances, not a proof.
- The formal core (operator content is RG-invariant up to symmetry-allowed mixing; symmetries are preserved by the
  flow) is standard. So the no-go reduces to "symmetry content is preserved by RG", which is **REDISCOVERED-KNOWN**
  in form and **NEW-IN-GRUT** as the boundary of the quotient program.

**Refined quotient principle (SCOUT-0 + SCOUT-1, candidate):**

> **Dynamics can erase parameter values but cannot choose its own operator content.** Observable IR quotients are
> functions of (operator / symmetry content, supplied) × (fixed-point data, universal) × (scales, G-moved, never
> fixed).

**Status: D1-SCOUT COMPLETE — B (operator class supplied) at the D1 level; C (KPZ forbidden) within inherited
structure; the pin reveals, does not select; no GRUT-internal selection. TC-4 SELECTOR RELOCATION NO-GO recorded as
CANDIDATE.**

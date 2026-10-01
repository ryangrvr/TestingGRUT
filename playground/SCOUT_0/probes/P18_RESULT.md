# SCOUT_0 W2 P-18 RESULT — the ergodic fast-variable null (one pass, no simulation)

**Charter:** `PROBE_CHARTERS.md` P-18 (frozen: fast–slow, quartic slow drift, no simulation;
success = a theorem-grade statement on whether martingale noise survives). It is framed by P-17:
**does a chaotic, non-harmonic deterministic bath escape P-17's class 𝓗 and restore
identifiability?** Script: `p18_homogenization_checks.py` (exact sympy; log
`p18_homogenization_checks.log`). NEW HYPOTHESIS CLASS — DOES NOT ALTER OLD TERMINAL.

## 0. Verdict

> **Under the declared weak-invariance/homogenization hypotheses, centred diffusively scaled forcing
> converges to Brownian noise whenever the asymptotic variance `Σ²` (Green–Kubo) is nonzero; no
> detailed-balance condition is needed. Coboundaries give `Σ² = 0` (exact example below); an "iff"
> characterization of `Σ² = 0` by coboundaries needs the relevant dynamical-class hypotheses (it holds
> in the classical hyperbolic settings) and is not claimed for arbitrary chaos. The limit law is C-B
> itself, with `Q = Σ²`. At every finite timescale separation the chaotic driver is, by P-17's
> argument, indistinguishable from an exogenous process with the same law. It does not escape
> non-identifiability. It is an AUTONOMOUS DETERMINISTIC CHAOTIC FORCING PARENT (a one-way driver),
> not a second physical bath.**
>
> *Repair note 01 (auditor, post-banking): three phrases repaired — the "exactly when … not a
> coboundary" statement (now hypothesis-scoped, above); "second physical parent/bath" (now
> "deterministic chaotic driver", §2–3); the Kelly–Melbourne citation (§1b).*

- **Charter hypothesis** ("not a true martingale unless detailed balance"): **FALSE.** The
  homogenization theorems need mixing plus a weak invariance principle, not detailed balance.
- **Charter null** ("no noise survives averaging"): **TRUE only in the averaging scaling**
  (`ẋ = f + h(y)`, fast `y`), or when `Σ² = 0` in the diffusive scaling (e.g. a coboundary `h`; exact
  example below). **FALSE** in the diffusive scaling whenever `Σ² > 0` under the WIP hypotheses.

## 1. Parent and statement

**Slow:** `ẋ = f(x) + ε⁻¹h(y)` with `f = −K_b x − 4βx³` (the C-B drift).
**Fast:** `ẏ = ε⁻²g(y)` (flow), or the discrete analogue `x_{n+1} = x_n + ε²f(x_n) + εh(y_n)`,
`y_{n+1} = T(y_n)`. The fast system is deterministic, mixing, with invariant measure `μ` (uniformly
or non-uniformly hyperbolic with summable correlations), `∫h dμ = 0`. Preparation: `x(0) = a·e₁`,
`y(0) ~ μ`.

**(a) Finite ε — exact non-identifiability.** There is no back-reaction (skew product). So the
forcing `F_ε(t) = ε⁻¹h(y(t/ε²))` has a law fixed by `y(0) ~ μ` and independent of `x(0)`. The slow
path is `Φ(x(0), F_ε)` for one solution map `Φ`. P-17's proof applies verbatim, with zero memory
kernel: **the reduced law equals that of `ẋ = f(x) + F` with `F` exogenous of law `Law(F_ε)`.**
Exact, at every ε, any `h`, Gaussian or not.

**(b) ε → 0 — C-B (IMPORTED-STANDARD).** Melbourne–Stuart (Nonlinearity 24, 2011),
Gottwald–Melbourne (Proc. R. Soc. A 469, 2013) and Kelly–Melbourne, "Deterministic homogenization
for fast–slow systems with chaotic noise", J. Funct. Anal. 272 (2017) 4063–4102 (the broad
deterministic-homogenization result; not the separate 2016 Ann. Probab. paper on smooth approximation
of SDEs), give, under their weak-invariance-principle hypotheses, weak convergence in path space:

`x_ε ⇒ X`, with `dX = f(X)dt + Σ dW` and `Σ² = C₀ + 2Σ_{n≥1}C_n`, `C_n = ∫h·h∘Tⁿ dμ`.

Noise is additive, so there is no Itô/Stratonovich ambiguity. For C-B, take one independent fast
system per site with `Σ_i² = 2T_i`.

**(c) Exact Green–Kubo checks (C1).** Doubling map `T(y) = 2y mod 1` (Lebesgue-invariant, mixing;
its natural extension, the baker's map, is invertible and area-preserving):

| `h` | `C₀, C₁, C₂, …` | `Σ²` | limit |
|---|---|---|---|
| `cos 2πy` | `1/2, 0, 0, …` | `1/2` | Brownian |
| `cos 2πy + cos 4πy` | `1, 1/2, 0, …` | `2` | Brownian |
| `χ∘T − χ`, `χ = cos 2πy` (coboundary) | `1, −1/2, 0, …` | `0` | **no noise survives** |

**(d) What reaches the C-B discriminator at finite ε (C2, C3).**

- The skew forcing `cos 2πy + cos 4πy` has `E h³ = 3/4`. The third cumulant of the slow increment is
  `O(ε)`, so the even-in-`a` term `−4β∫g κ₃` (P-17 lemma) vanishes in the limit. The odd-in-`a` part
  is covariance-only and tends to the C-B value with `T₁ = Σ₁²/2`.
- **Regime fragility (same as P-17):** for smooth forcing, at `t ≪ ε²τ_c` the leading term is
  `S − D ≈ −4βaσ_F²t³`, not C-B's `t²`. The C-B coefficient holds in the window `ε²τ_c ≪ t`.

## 2. Ledgers

| Ledger | Answer |
|---|---|
| PHYSICAL PARENT EXISTS | **An autonomous deterministic chaotic forcing parent exists — not a physical bath.** It is a one-way driver of C-B (non-Gaussian at finite ε), at homogenization-limit grade. The skew product has **no back-reaction**, so the combined slow–fast system is not a reciprocal, energy-conserving Hamiltonian bath, even when the fast subsystem alone is Hamiltonian or area-preserving (Anosov geodesic flow; baker's map). Contrast: P-17 = reciprocal Hamiltonian environment; P-18 = one-way deterministic chaotic driver. |
| REDUCED PROCESS LOOKS STOCHASTIC | Under `y(0) ~ μ`: yes. For a single `y(0)`: deterministic path. In the limit, Brownian. |
| TRUE INNOVATIONS IDENTIFIABLE | **No** — (a) at finite ε, (b) in the limit. Escaping would require **back-reaction**, where the bath's state depends on the system path beyond a deterministic memory functional. Even there, the known homogenization limits (state-dependent drift and diffusion corrections) are again Markov diffusions that an exogenous multiplicative-noise model reproduces. |

## 3. Structural reading (feeds ZOOM-OUT 02)

- **Universality.** Two unrelated microscopic constructions give the same reduced law, **C-B**:
  1. P-17: a **reciprocal Hamiltonian environment** (Gaussian harmonic bath, FDT), through the WB + OD
     limits;
  2. P-18: a **one-way deterministic chaotic driver** (non-Gaussian at finite ε), through
     homogenization.

  Only one nuisance parameter survives per site: `Q = 2T` (FDT temperature in one parent, Green–Kubo
  integral in the other). This is the standard CLT/invariance-principle universality
  (**REDISCOVERED-KNOWN**). In GRUT terms: **the C-B class is the universal diffusive limit, and its
  discriminator's leading odd part reads only `T₁`.**
- **Non-identifiability is generic, not bath-specific.** P-17 and P-18 share one mechanism: whenever
  the force on the retained system splits into (a force fixed by hidden initial data, with a law
  independent of the system preparation) plus (a deterministic functional of the system path), the
  reduced law equals that of an exogenous model. With S2-HB's HB-U path-space realization (which
  always exists), **no function of a reduced law can certify primitive randomness**. Reduced data
  identify the *law*; which ontology realizes it is a question about admissible physical classes,
  not about observations.

**Classification:** REDISCOVERED-KNOWN (homogenization; Green–Kubo; coboundary degeneracy).
KNOWN-BUT-NEW-IN-GRUT: a deterministic chaotic driver of C-B (contrasted with P-17's reciprocal
Hamiltonian environment), and the generic form of P-17's non-identifiability.

**Status: P-18 COMPLETE (one pass; repair note 01 applied). Under WIP/homogenization hypotheses,
centred diffusive forcing with `Σ² > 0` gives Brownian noise without detailed balance; coboundaries
give `Σ² = 0`; limit = C-B with Green–Kubo `Q`; exact non-identifiability at every ε; an autonomous
deterministic chaotic driver of C-B (not a bath). Old terminals untouched.**

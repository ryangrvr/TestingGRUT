# SCOUT_0 W2 P-18 RESULT — the ergodic fast-variable null (one pass, no simulation)

**Charter:** `PROBE_CHARTERS.md` P-18 (frozen: fast–slow, quartic slow drift, no simulation;
success = a theorem-grade statement on whether martingale noise survives). It is framed by P-17:
**does a chaotic, non-harmonic deterministic bath escape P-17's class 𝓗 and restore
identifiability?** Script: `p18_homogenization_checks.py` (exact sympy; log
`p18_homogenization_checks.log`). NEW HYPOTHESIS CLASS — DOES NOT ALTER OLD TERMINAL.

## 0. Verdict

> **Noise survives as a true martingale (Brownian) limit without any detailed-balance condition,
> exactly when the forcing is centred, diffusively scaled, and not a coboundary. The limit law is
> C-B itself, with `Q = Σ²` given by Green–Kubo. At every finite timescale separation the chaotic bath
> is, by P-17's argument, indistinguishable from an exogenous process with the same law. It does not
> escape non-identifiability; it gives a second, independent physical parent of C-B.**

- **Charter hypothesis** ("not a true martingale unless detailed balance"): **FALSE.** The
  homogenization theorems need mixing plus a weak invariance principle, not detailed balance.
- **Charter null** ("no noise survives averaging"): **TRUE only in the averaging scaling**
  (`ẋ = f + h(y)`, fast `y`), or for a coboundary `h` in the diffusive scaling (exact example below).
  **FALSE** in the diffusive scaling for generic `h`.

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
Gottwald–Melbourne (Proc. R. Soc. A 469, 2013) and Kelly–Melbourne (Ann. Probab. 44, 2016) give
weak convergence in path space:

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
| PHYSICAL PARENT EXISTS | **Yes — a second, structurally different parent of C-B** (deterministic chaos, non-Gaussian at finite ε), at homogenization-limit grade. It is energy-conserving only when the fast flow is Hamiltonian chaos (e.g. Anosov geodesic flow, or the baker's map in discrete time); the skew product has **no back-reaction**, which fails a strict "physical bath" reading (action–reaction). That is the price, recorded. |
| REDUCED PROCESS LOOKS STOCHASTIC | Under `y(0) ~ μ`: yes. For a single `y(0)`: deterministic path. In the limit, Brownian. |
| TRUE INNOVATIONS IDENTIFIABLE | **No** — (a) at finite ε, (b) in the limit. Escaping would require **back-reaction**, where the bath's state depends on the system path beyond a deterministic memory functional. Even there, the known homogenization limits (state-dependent drift and diffusion corrections) are again Markov diffusions that an exogenous multiplicative-noise model reproduces. |

## 3. Structural reading (feeds ZOOM-OUT 02)

- **Universality.** Two unrelated microscopic parents give the same reduced law, **C-B**:
  1. P-17: a Gaussian harmonic Hamiltonian bath, with FDT, through the WB + OD limits;
  2. P-18: non-Gaussian deterministic chaos, through homogenization.

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
KNOWN-BUT-NEW-IN-GRUT: the second physical parent of C-B, and the generic form of P-17's
non-identifiability.

**Status: P-18 COMPLETE (one pass). Martingale noise survives without detailed balance (non-coboundary,
diffusive scaling); limit = C-B with Green–Kubo `Q`; exact non-identifiability at every ε; a second
universal parent of C-B. Old terminals untouched.**

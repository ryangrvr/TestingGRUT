# V0-1 TARGET SPEC — P-17 (statement, assumptions, definitions, acceptance tests ONLY; no proof)

**Source:** frozen `scout-0 @ ab2da47`:
- `playground/SCOUT_0/probes/P17_RESULT.md`, blob `5fc0401b60de2e91f787acb7fa1ae4c55d70f396`;
- the acceptance-test parameters and the preparation definitions come from `p17_taylor_identity.py` (lines 12 – 56 read for
  parameters / preparations only) and from its logs `p17_taylor_2modes_n8.log` and `p17_taylor_1mode_n12.log`.

**After this spec is committed, those files are SEALED until the V0-1 reproduction is committed.**

## 1. Class 𝓗 (definitions)

**System.** q ∈ ℝⁿ, momentum p, H_S = ½pᵀM⁻¹p + V(q), with V smooth and confining (any nonlinearity).

**Bath and coupling.** Harmonic modes (x_j, p_j), finitely many. The coupling is linear in x_j and in q, with a
counterterm:

  H = H_S + Σ_j [ p_j²/2 + ½ ω_j² (x_j − c_jᵀq/ω_j²)² ]   (scalar modes; c_j ∈ ℝⁿ).

**Preparation hypothesis (P).** Given the system initial state (q₀, p₀), the free-force process F (the part of the bath
force determined by the bath initial data alone; its precise form is to be derived) has a law 𝒫 that does **not**
depend on (q₀, p₀).

## 2. Target claims (to be re-derived independently)

**THEOREM (P-17).** Let A be the Hamiltonian system with bath initial states drawn per (P). Let B be the GLE

  M q̈ = −∇V(q) − ∫₀ᵗ γ(t − s) q̇(s) ds + F(t),  with  γ(t) = Σ_j (c_j c_jᵀ/ω_j²) cos ω_j t,

driven by an **exogenous** process F ~ 𝒫, independent of (q₀, p₀). Then **for every (q₀, p₀), the laws of the system
path q(·) on any [0, t_max] coincide under A and B.** Hence every reduced statistic agrees.

**Claimed consequences:**
- **(i) Thermal shifted-Gibbs preparation** (bath Gibbs at temperature T around x_j = c_jᵀq₀/ω_j²): F is Gaussian with
  mean 0 and covariance T·γ(t − s). Matched covariance suffices.
- **(ii) Product preparation** (bath Gibbs around x_j = 0, ignoring q₀): (P) fails. A equals B plus a deterministic
  initial-slip forcing **−γ(t)q₀**, exactly.
- **(iii) Deterministic bath microstate:** 𝒫 is a Dirac law, and the identity is pathwise.
- **(iv) Ontology conclusion:** reduced data identify the effective forcing **law** (its cumulant hierarchy) and the
  memory kernel. They never identify whether the forcing came from hidden deterministic initial conditions or from
  primitive randomness, **within class 𝓗**.
- **(v) C-B interpretation:** the record's C-B discriminator reads the forcing-law cumulants. It detects forcing from
  outside the retained state space, not ongoing randomness. The exact C-B (white, overdamped) form is reproduced only in
  the non-admitted white-bath + overdamped limits.

## 3. Acceptance tests (exact rational Taylor checks; "C2")

**Model.**
- n = 1, M = 1, V(q) = (23/10)·q²/2 + β q⁴.
- Modes (ω_j, c_j): mode 1 = (1, 1/2), mode 2 = (2, 3/4).
- **Two-mode** runs use both modes, through Taylor order n ≤ 8. **One-mode** runs use mode 1 only, through n ≤ 12.

**Initial data.** q(0) = a, p(0) = 0, with x_j(0) = c_j a/ω_j² + ξ_j and p_j(0) = η_j.

**Preparations:**
- **Thermal (A):** ξ_j ~ N(0, T/ω_j²) and η_j ~ N(0, T), all independent.
- **Random-phase (R):** ξ_j = (√(2T)/ω_j) cos φ_j and η_j = −√(2T) sin φ_j, with φ_j i.i.d. uniform.
- **Product:** bath Gibbs around x_j = 0 regardless of q₀, i.e. x_j(0) ~ N(0, T/ω_j²) and p_j(0) ~ N(0, T).

**Comparators:**
- **B:** the GLE above with exogenous Gaussian F, mean 0, covariance T·γ.
- **D:** the deterministic GLE, F ≡ 0.
- **B_R:** the GLE with exogenous F having the same law as the random-phase free force (independent phases).

**Statistic:** E[q⁽ⁿ⁾(0)], exact as polynomials in (β, T, a).

| test | expected result |
|---|---|
| **T1** (A vs B) | A − B = 0 for all n ≤ 8 (two modes) and n ≤ 12 (one mode) |
| **T2** (A − D), two modes | 0 for n ≤ 5 and all odd n. n = 6: −225·T·a·β/8. n = 8: 15·T·a·β·(63360·a²β + 25973)/256 |
| **T2** (A − D), one mode | 0 for n ≤ 5 and all odd n. n = 6: −18·T·a·β. n = 8: 24·T·a·β·(495·a²β + 193)/5. n = 10: −3·T·a·β·(27302400·a⁴β² + 18202640·a²β + 1976871)/200. n = 12: 3·T·a·β·(46170000·Tβ + 15956352000·a⁶β³ + 15791424000·a⁴β² + 4083013860·a²β + 136763237)/500 |
| **T3** | the leading A − D coefficient equals 6!·[−β·a·T·γ(0)/(10M³)] |
| **T4** (A vs R) | A − R = 0 for all n ≤ 8 (two modes). One mode: 0 for n ≤ 11, and at n = 12, A − R = 138510·T²·a·β² |
| **T5** (R vs B_R) | R − B_R = 0 at every computed order (both configurations) |
| **T6** (B_R vs B) | zero for n ≤ 8 (two modes). One mode: non-zero only at n = 12 |
| **T7** (product) | product − B (shift-matched) is non-zero at even n: n ∈ {2, 4, 6, 8} for two modes, {2, …, 12} for one mode. Product − (B + deterministic slip −γ(t)a) = 0 at all computed n |
| **T8** (controls) | β = 0 ⇒ A − D = 0 at all n. T = 0 ⇒ A − D = 0 at all n |

## 4. Scope / exclusions

- **Exactness needs (P)** (or a violation absorbable into deterministic memory, as in (ii)).
- **Outside 𝓗** (anharmonic bath modes, or coupling nonlinear in the bath coordinates), the theorem is **not claimed**.
- **Classical only.** The quantum noise-kernel remarks are out of V0-1 scope.
- **Not in V0-1:** the C3 short-time regime table and the C4 white-C-B cumulant identities (Δc₂ = −24βTa − 4βκ₃; κ₄
  entering Δc₃ as +240β²κ₄a). They need the full C-B model definition from elsewhere in the record. Registered as
  deferred, not owed by V0-1.

**Hostile scope checks required by the owner:**
- non-harmonic bath;
- coupling nonlinear in q;
- a q₀-dependent preparation;
- extra restrictions on the exogenous force.

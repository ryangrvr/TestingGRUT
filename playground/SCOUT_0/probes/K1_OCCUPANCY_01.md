> **CORRECTED — see `K1_OCCUPANCY_CORRECTION_01.md` and `K1H_RESULT.md`.** Two conclusions below are WITHDRAWN:
> (1) "Only Run sits exactly on Σ₀ = μ₀/2 at z = 0 for every amplitude" — Only Run gives `2Σ − μ = R`, so it is on K1
> only where `M_*² = m_p²`, which in Linder's benchmark is the early GR-restoration epoch, not z = 0;
> (2) "holding K1 at all a forces GR within leading-order Horndeski" — the exact condition is a first-order ODE that
> admits non-GR solutions (K1-H: stable ones exist on the μ < 1 half-line). Linder 2020 formulas upgraded to
> PRIMARY-TEXT-VERIFIED (external audit). The text below is kept unchanged for the audit trail.

# SCOUT_0 K1 OCCUPANCY CHECK 01 — does an established theory occupy `Σ₀ = μ₀/2`?

**Purpose:** a hostile discriminator check on K1 (`P24_RESULT.md`) before anyone considers the fenced
`p_tt`/Π₀ derivation. **Question:** would nature distinguish K1 from established
modified-gravity/dark-energy families if GRUT derived it? Read-only plus literature. Script for the
algebra: `k1_occupancy_checks.py` (log `k1_occupancy_checks.log`).

**Source grade:** arXiv, A&A, OSTI and INSPIRE are blocked by this sandbox's egress proxy. Literature
relations below are therefore **IMPORTED-SECONDARY** (search-engine excerpts of the primary papers).
They are cross-checked against independently known limits where possible (marked ✓). Each owes a
primary-text check before any external use.

## 1. GRUT conventions, reconstructed from the record

| Item | Record definition | Source |
|---|---|---|
| `μ` | growth: Poisson modifier of the Newtonian potential Ψ (`μ = 1 + xα`, "exact in the inherited bookkeeping") | `GRUT_ToE.md` §2.8; `GRUT_II_What_Survived.md` §5 |
| `η` | slip `Φ/Ψ`; family value `η = 1/(1 + xα)` | same; `calc/sigma0_anomaly_screen.py` |
| `Σ` | lensing, `Σ = μ(1+η)/2`; family value `1 + xα/2` | same |
| `μ₀, Σ₀` | GRUT: `μ − 1`, `Σ − 1`, **constant in time and scale** on the constant-x cut. The observational papers use `μ(a) = 1 + μ₀Ω_DE(a)/Ω_Λ` (Ω_DE-tracking), so amplitudes compare directly only at `z = 0` | `calc/sigma0_anomaly_screen.py` ("F-MAP fence") |
| **What the lock is** | `Σ − 1 = (μ − 1)/2` ⟺ **`μη = 1`** ⟺ the curvature potential Φ is unmodified, only Ψ is modified. It holds pointwise wherever `η = 1/μ` holds, so it is a relation between functions of `(k, a)`, not only of amplitudes | `k1_occupancy_checks.log` |
| Conditionality | `η = 1/μ` (hence the lock) is **INHERITED-CONDITIONAL** on the bookkeeping's one-sidedness: "a genuine P⁰ˢ admixture could split the law via scalar anisotropic stress (R2, open)". Only the μ law is unconditional. The constant-x family is "a one-dimensional cut through a function space" (x is a kernel) | `GRUT_II_What_Survived.md` §5; `S_IF.md` |

So K1 has two readings: **K1-amp** (the `z = 0` amplitude relation `Σ₀ = μ₀/2`) and **K1-fun** (`μη ≡ 1`
at all `(k, a)`, giving a time-constant ratio `(Σ−1)/(μ−1) = 1/2`).

## 2. Comparison classes

| Class | μ–Σ relation (quasi-static) | Time / scale dependence | Holds when | On GRUT line? |
|---|---|---|---|---|
| ΛCDM / GR | `μ = Σ = η = 1` | — | always | origin only (GRUT's `x = 0`) |
| f(R) (e.g. Hu–Sawicki) | `Σ = 1/F ≈ 1`; `μ = (1+4κ)/(1+3κ) ∈ [1, 4/3]` | scale-dependent (Compton) | quasi-static | **no** (`Σ₀ ≈ 0`), except origin ✓ (record: `calc/anomaly_c0_map.py`) |
| Brans–Dicke / conformal scalar-tensor | `μ = (2ω+4)/(2ω+3)`, `η = (ω+1)/(ω+2)` ⇒ `Σ = 1`; `μη = 2(ω+1)/(2ω+3) ≠ 1` | constant in ω | light field, quasi-static | **no**, except origin (algebra ✓) |
| DGP (normal branch) | brane-bending mode conformally coupled ⇒ `Σ = 1`, `μ = 1 + 1/(3β)` | time-dependent β | quasi-static | **no**, except origin |
| No Slip Gravity (Linder 2018) | `η = 1` ⇒ `Σ = μ` | free α-functions | Horndeski subclass, `α_T = 0` | **no** (line `Σ₀ = μ₀`), except origin |
| No Run Gravity / cubic Galileon / KGB (`α_M = 0`) | no slip ⇒ `Σ = μ` | — | cubic Horndeski | **no**, except origin |
| General Horndeski, `α_T = 0` (post-GW170817), small α, `M_* ≈ m_p` | `(Σ−1)/(μ−1) ≈ (1 + r)/(2 + r)`, `r = α_B/α_M`. Limits ✓: f(R) `r = −1` → 0 (Σ = 1); `α_M → 0` → 1 (no slip) | functions of `a` | leading order in α | **only on the locus `α_B = 0`**. Pogosian–Silvestri conjecture `(Σ−1)(μ−1) ≥ 0` (confirmed by Peirone et al.): the GRUT line lies **inside** the Horndeski-allowed region, so it is not excluded |
| **Only Run Gravity** (Linder 2020, `α_B = 0`) | `G_matter = R(1+α_M)`, `G_light = R(1+α_M/2)`, `R = m_p²/M_*²` ⇒ **`Σ−1 − (μ−1)/2 = ε/2`**, `ε = R − 1`; `μη = R` | `α_M(a)` free; `ε(a)` = accumulated Planck-mass running | quasi-static, `α_T = 0` | **ON the line wherever `M_* = m_p`** (the normalization epoch), for **every** `α_M` value there — a full one-parameter intersection. Off it elsewhere by `ε(a)/2`. Staying on it at all `a` forces `ε ≡ 0 ⇒ α_M ≡ 0` ⇒ GR |
| Phenomenological μ₀–Σ₀ (DES, KiDS, Du et al.) | independent amplitudes | Ω_DE-tracking | parameterization, not a theory | the line is a 1D subset; no occupant |
| Coupled quintessence (DM-only coupling) | metric potentials GR-like (`Φ = Ψ`, `Σ = 1`); the fifth force shows up as a growth `μ_eff = 1 + 2β²` | — | quasi-static | **no** (`Σ₀ = 0` in lensing), except origin (convention caveat: growth-only `μ_eff` is not a metric μ) |
| Clustering DE (`c_s → 0`, no anisotropic stress) | `η = 1` ⇒ `Σ = μ` | — | — | **no**, except origin |
| Beyond-Horndeski / DHOST, nonlocal (Deser–Woodard, Maggiore RT/RR), massive/bigravity, DE with anisotropic stress | **NOT EVALUATED** in this pass (more free functions; may reach the line by tuning) | — | — | open |

## 3. Verdict

> **K1-amp (`z = 0` amplitudes): C. POINT-INTERSECTION — but a dangerous one.**
> - Only Run Gravity (`α_B = 0`, running Planck mass, no braiding) sits **exactly** on
>   `Σ₀ = μ₀/2` at its normalization epoch for every amplitude.
> - Low-redshift amplitude-only data therefore **cannot** distinguish GRUT's lock from Only Run.
> - Every other compared class meets the line only at the origin (`Σ₀ = 0` or `Σ₀ = μ₀`).
>
> **K1-fun (`μη ≡ 1` at all `(k, a)`; time-constant ratio 1/2): A. DISTINCT-LINE within the declared
> comparison set.**
> - Only Run departs from it by `ε(a)/2`.
> - Within leading-order `α_T = 0` Horndeski, holding the line at all `a` forces GR.
> - The non-evaluated classes (beyond-Horndeski, DHOST, nonlocal, massive gravity) remain open.
> - **"No occupancy found in the declared comparison set"** — not a uniqueness claim.

**Consequence for the `p_tt` decision.** Deriving `p_tt` would **not** by itself produce a distinctive
prediction at the level of today's amplitudes; Only Run Gravity shares that. The distinctive content
is the **time dependence**: GRUT's lock holds at every redshift (Φ unmodified throughout), while the
nearest established occupant drifts off it by half the accumulated Planck-mass running.

So a `p_tt` derivation is worth pursuing **only together with** the redshift-resolved test:
`μη(z) = 1` (GRUT) vs `μη(z) = m_p²/M_*²(z)` (Only Run). Two further conditions apply:
- the derivation must also settle the inherited-conditional slip `η = 1/μ` (R2), without which
  K1-fun itself is not GRUT content;
- the observational separation is a second-order effect (`ε` is an integral of `α_M`), so its
  reachable precision is an owed estimate.

**Owed:**
- (i) primary-text verification of the Only Run and small-α Horndeski formulas;
- (ii) evaluation of the beyond-Horndeski, DHOST and nonlocal classes;
- (iii) a forecast-level estimate of the `μη(z)` separation, GRUT vs Only Run, at plausible survey
  precision.

## Sources (literature excerpts used)

- Pogosian & Silvestri, "What can Cosmology tell us about Gravity? Constraining Horndeski with Sigma and
  Mu", arXiv:1606.05339 — https://arxiv.org/abs/1606.05339 (conjecture `(Σ−1)(μ−1) ≥ 0`)
- Peirone et al. (confirmation of the conjecture) — via the excerpt at
  https://arxiv.org/pdf/1902.10503
- Linder, "No Slip Gravity", arXiv:1801.01503 — https://arxiv.org/abs/1801.01503
- Linder, "No Run Gravity", arXiv:1903.02010 — https://arxiv.org/abs/1903.02010
- Linder, "Limited Modified Gravity", JCAP 10 (2020) 042, arXiv:2003.10453 —
  https://arxiv.org/abs/2003.10453 (Only Run: `G_matter`, `G_light`)
- Small-α Horndeski ratio `(1 + α_B/α_M)/(2 + α_B/α_M)` — excerpts from https://arxiv.org/pdf/2011.01517
  ("Weak gravity on a ΛCDM background") and KiDS-Legacy Horndeski constraints
  (https://arxiv.org/pdf/2512.11039)

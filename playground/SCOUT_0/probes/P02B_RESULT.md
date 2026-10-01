# SCOUT_0 W3 P-02b RESULT — interacting bosons: does density → IR velocity replace exclusion?

**Purpose:** attack the quadratic/free scope of the EDGE-DATA AUDIT theorem (`EDGE_DATA_AUDIT_01.md`
§4). P-02 (free bosons, BOSON-COLLAPSE) is the unmodified control.
**Parent:** `H = −Σ_j (b_j†b_{j+1} + h.c.) + (U/2)Σ_j n_j(n_j − 1)` on SF-1's periodic ring. The
interaction `U` is a **NEW ASSUMPTION**. Same conserved `N`, readout `ρ_q`, support object and
sector rules as SF-1. U = 0 is exactly P-02; U = ∞ is the hard-core diagnostic.
**Script:** `p02b_bose_hubbard.py` (log `p02b_bose_hubbard.log`). The sector ground state is unique
for every U (Perron–Frobenius: hopping off-diagonals ≤ 0, connected configuration graph; checked
gaps > 0). The lowest supported excitation at `q` = lowest Ritz value of `H` on the Krylov space of
`ρ_q|0⟩` (full reorthogonalization). Exact sector Hilbert spaces up to L = 14 (dimension 77 520).

## 0. Verdict

> **Grade note (external audit):** the statement "for every U > 0" is **IMPORTED-STANDARD** 1D
> Bose–Hubbard / Luttinger / Tonks physics. The L ≤ 14 exact diagonalization is **finite-size
> confirmation**, not the proof, and is not presented as one anywhere below.
>
> **SAME-SPLIT restored by interaction (1D).** For every U > 0 the dense sector (ν = 1/2) has `z → 1`
> and a soft point at `q = 2πν = π`: class `(1, 2)`, as for SF-1's fermions. The dilute sector
> (fixed N) has `z → 2`: class `(2, 1)`.
>
> - **Exclusion is not necessary:** a repulsive interaction supplies the density → IR velocity map
>   instead.
> - **The free-boson collapse (P-02) is a singular point at U = 0**, not a generic bosonic property.
>   The class jumps discontinuously at U = 0⁺.
>
> **Consequence for the audit theorem:**
> - **The "values are supplied" half survives** (U and ν are supplied).
> - **The "factors through one-particle edge data" half FAILS** for interacting parents. The
>   relevant velocity is a many-body (Luttinger/sound) velocity `v(ν, U)`, not `ε′(k_b)`.
> - The soft-point count additionally depends on the **IR fixed-point type** (1D Luttinger: 2 soft
>   points; condensate/Bogoliubov: 1). That is genuinely new law-class data beyond the free edge
>   dictionary.

## 1. Results

**D-type sector** (`N = L/2`, L = 6, 10, 14; `ω₁ = ω⁻(2π/L)`; local `z_{P′}` between successive L):

| U | `ω₁(L)` | local `z_{P′}` | `ω⁻(q = π)` | reading |
|---|---|---|---|---|
| 0 (control) | 1.000, 0.382, 0.198 | 1.88, 1.95 → 2 | 4.000, 4.000, 4.000 | collapse: band bottom, no soft point at π |
| 1 | 1.282, 0.660, 0.444 | 1.30, 1.18 → 1 (slow: healing-length crossover) | 2.92, 1.87, 1.36 (→ 0) | Luttinger |
| 4 | 1.717, 1.013, 0.719 | 1.03, 1.02 | 2.62, 1.69, 1.23 (≈ 1/L) | Luttinger |
| ∞ (hard-core) | 2.000, 1.236, 0.890 | 0.94, 0.98 → 1 (exact: = SF-1 D) | 2.000, 1.236, 0.890 | free-fermion map |

- Full edge at L = 10, U = 4: `ω⁻(q) = 1.01, 1.69, 2.03, 2.04, 1.69, …` (symmetric). There is a local
  dip at `q = π` (m = 5) between maxima, the finite-L shadow of the 2k_F soft point. It sharpens with
  L (column above).

**E-type sector** (fixed `N = 3`, L = 6, 10, 14):

| U | `ω₁(L)` | local `z_{P′}` |
|---|---|---|
| 0 | 1.000, 0.382, 0.198 | 1.88, 1.95 → 2 |
| 4 | 1.717, 0.758, 0.429 | 1.60, 1.69 (rising toward 2) |
| ∞ | 2.000, 1.000, 0.555 | 1.36, 1.75 → 2 (exact: `ω₁ = 4 sin(π/L) sin(3π/L)`, SF-1 E-3) |

**Dilute-limit argument (why E stays z = 2 for any U).** At fixed N and L → ∞ the density
`n → 0`. The Lieb–Liniger coupling `γ = U/n → ∞` drives the gas to the Tonks (hard-core) limit,
which is exactly free fermions with `k_F → 0`: a band-bottom edge, `z = 2`. The U = 4 numerics are
still on that crossover at L ≤ 14 (1.60 → 1.69). **The E-sector z = 2 for finite U rests on this
standard limit argument; the L ≤ 14 numerics show the trend, not the limit.**

**Mean-field / d ≥ 2 contrast (analytic Bogoliubov).** `ω_B = √(ε_q(ε_q + 2Un))`:
- finite n gives `ω ≈ √(2Un) q` (z = 1) with soft set `{0}` only → class **(1, 1)**;
- n → 0 gives `ω → ε_q` (z = 2) → (2, 1).

So a condensate splits z the same way but has **one** soft point, while 1D Luttinger liquids have two
(the Lieb II mode reaches 0 at `2πn`). Bogoliubov is not the exact 1D answer; it stands in for the
condensate IR type.

## 2. Minimal carrier, updated

| Mechanism | Density → IR velocity? | D vs E split in z? | Soft points (D) |
|---|---|---|---|
| exclusion (fermions, hard-core bosons) | `v = ε′(πν)` (one-body edge) | yes | 2 |
| repulsive interaction, 1D | `v = v_s(ν, U)` (many-body, Luttinger) | yes | 2 |
| repulsive interaction, condensate (mean-field / d ≥ 2) | `v = √(2Un)` (Bogoliubov) | yes | 1 |
| none (free bosons, U = 0) | `v ≡ 0` | **no** (collapse) | 1 |

> **Refined statement (candidate):** within these parents, the sector → law-class map factors through
> **(i)** the IR velocity of the sector's excitation edge, `v_IR(ν; supplied mechanism)` (z = 1 iff
> `v_IR ≠ 0`), and **(ii)** the IR fixed-point type (Luttinger vs condensate), which sets the
> soft-point count. Conservation is needed but not sufficient. Exclusion and interaction are
> interchangeable suppliers of (i). Neither is earned.

## 3. What this does to the audit theorem and the conjecture

- **Audit theorem** (`EDGE_DATA_AUDIT_01.md` §4) stays valid **at its declared quadratic/free
  scope**. P-02b shows the scope matters: with interactions, the quotient is not one-particle edge
  data. It is **many-body IR data** (an excitation-edge velocity plus an IR fixed-point type).
- **Supplied-values half: survives.** `U`, ν and dimension are supplied; nothing earned fixes `v_IR`
  or the IR type.
- **Conjecture "GRUT restricts the law iff it fixes edge data":** survives only if "edge data" means
  **many-body IR data**, including the fixed-point type. In that form it is close to a restatement of
  RG universality (IR class = fixed point + relevant data). **Recorded as REFINED, still a
  conjecture.** It is not promoted.
- **Non-genericity of P-02's collapse:** the free-boson collapse holds only at U = 0 (the class is
  discontinuous there). The P-02 result is not modified; its reading gains this qualifier.

**Classification:** REDISCOVERED-KNOWN (Luttinger-liquid sound velocity; Lieb–Liniger Tonks limit;
Bogoliubov phonons; Lieb II 2k_F softening). KNOWN-BUT-NEW-IN-GRUT: interaction replaces exclusion as
the SF-1 carrier, the free-boson collapse is singular at U = 0, and the edge quotient must be lifted
from one-body spectral data to many-body IR data.

**Status: P-02b COMPLETE (one controlled pass). SAME-SPLIT for U > 0 in 1D (IMPORTED-STANDARD Luttinger/Tonks physics; L ≤ 14 ED = finite-size confirmation only) (D: z → 1 with a soft point
at π; E: z → 2, via the Tonks argument for finite U). The free-boson collapse is singular at U = 0.
The audit theorem holds at its quadratic scope only; the interacting quotient is many-body IR data
(velocity + fixed-point type), still supplied. Next in sequence: P-15.**

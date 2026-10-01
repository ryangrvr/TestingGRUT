# SCOUT_0 W3 P-24 RESULT — lock census (read-only)

**Purpose (ZOOM_OUT_02):** attack P-23's enumeration-grade claim that `Σ₀ = μ₀/2` is the record's only
surviving law-level lock. A **lock** = a relation between two observables in which every supplied
amplitude/parameter cancels. Each lock found is classified: does it survive variation of the
supplied inputs, does it bear on an empirical observable, does it differ from a standard class?

**Method (search grade):**
- Pattern search over the root-level verdicts, owner rulings, ledgers and synthesis files, plus
  `provenance/claims.json`. Terms: `lock`, `cancels`, `α-free`, `amplitude-free`, `parameter-free`,
  `independent of …`, `ratio`, `exponent`, `ω⁷`.
- Excluded: `archive/`, `release/`, `uploads/`, `books/` (superseded or derivative).
- Every hit read in context. This is a census of what the record **banks**, not a derivation of new
  locks.

## Census

| # | Lock (record source) | Cancels | Empirical observable? | vs standard | Class |
|---|---|---|---|---|---|
| K1 | **`Σ₀ = μ₀/2`, `μη = 1`, time-constant** (`calc/sigma0_anomaly_screen.py`: "EXACT LOCK, inherited-conditional") | `x`, `α` | **yes** (growth vs lensing, `μ₀–Σ₀` plane) | **off** the conformal/f(R)/DGP line `Σ₀ = 0` | **CONDITIONAL LAW-LEVEL LOCK** (rests on CHOSEN `p_tt`; vacuous at `x = 0`) |
| K2 | Noise/dissipation `N/Im K_R = 2coth(ω/2T)`, `c₀ = α` cancels (`NO_GO_LEDGER.md` §1; claims `ledger_note`: "N locked to Im χ") | `c₀ = α` | in principle | **standard** (KMS/Callen–Welton; BORROWED) | STANDARD |
| K3 | GR constraint lock `P⁰ˢ/P² = −2` (`GRUT_II_What_Survived.md:36`; dispatch) | normalization | via the scalar/tensor kernel | **standard GR** (linearized Einstein operator) | STANDARD |
| K4 | `Γ_T` horn (a): coefficient ratio `104/9` and chromatic `ω³` (`GRUT_PREDICTION_GATE_GAMMA_T.md` §2) | amplitude-free given the μ-slot | GW friction | shape = 4D stress-tensor spectral weight `∝ ω⁴` (dimensional) | **INVISIBLE** (≳60 orders, R5) |
| K5 | `w` one-signed, no phantom crossing (`NO_GO_LEDGER.md` §3; prediction map §0) | amplitude-free (sign) | **yes** (DESI DR3 exposure) | shared with single-field quintessence and ΛCDM | forbidden region, **NON-DISTINCTIVE**; `to-derive` |
| K6 | Matsubara comb at spacing exactly `H` (prediction map P2) | amplitude-free | not established | spacing `2πT_dS = H` is standard dS (KMS/quasinormal) | **FENCED**, unrun; spacing standard even if inherited |
| K7 | α-family exports `R = √(1+α)`, `S = 12π/α²`, `Ω_Λ(α)` (claims `rung9a_value`) | eliminating `α` gives e.g. `S = 12π/(R²−1)²` | **not identified** (R is a FRONTIER-RESERVED founding hypothesis; S and Ω_Λ formulas not located as observables in this pass) | — | **UNIDENTIFIED-OBSERVABLE**; conditional on adopted α |
| K8 | `ω⁷` exponent class (GR1/GR2A/CP-1/EQ-1) | coupling normalization | toy-internal | "the core does not select ω⁷"; "exponent kinematic" | NON-EMPIRICAL, CONDITIONAL |
| K9 | Toy-instrument ratios: FS-1 `Var(T=0.5)/Var(T=2)` (amplitude-free, separates H/Q₁/P₁); G-2 speed/stiffness `1/2.25`; GR2-d2 velocity ratio `1.809` | amplitudes | within declared toy models | access/geometry diagnostics | NON-EMPIRICAL (internal discriminators) |
| K10 | Open-system locks from this campaign: P-23a variance relation; `Γ(m₁)/Γ(m₂)`; C-B odd part reads `T₁` only | bath amplitude/width | yes | **standard** (P-23a; DP/CSL) | STANDARD |

## Verdict

> **Census-confirmed (search grade): K1 `Σ₀ = μ₀/2` is the unique record lock satisfying the census
> criteria (law-level, parameter-eliminated, tied to an identified empirical observable) and distinct
> from the standard classes compared so far** — it is off the compared conformal/f(R)/DGP line
> `Σ₀ = 0`. Whether another established modified-gravity class occupies the same time-constant
> `Σ₀ = μ₀/2` line is an **owed literature question**; until it is answered K1 is not called "off
> a standard class" without that qualification. *(Wording repair 01, external audit.)*
>
> - Every other lock is standard (K2, K3, K10), invisible (K4), non-distinctive (K5), fenced (K6),
>   toy-internal (K8, K9), or lacks an identified observable (K7).
> - P-23's claim is upgraded from enumeration grade to **census grade**.

**Residual risk:**
- **K7** is the one item not fully classified. If `S` or `Ω_Λ(α)` is an identified observable with a
  standard comparison, an α-eliminated lock between them could join K1, still conditional on the
  adopted α.
- The search used fixed terms. A lock stated without any of those words would be missed.

**Owed follow-ups (owner's call, not opened):**
- (i) Locate the `S` and `Ω_Λ(α)` definitions and classify K7.
- (ii) A literature check of whether any standard theory occupies the line `Σ₀ = μ₀/2` with a
  time-constant shape. This decides whether K1, if ever derived, would discriminate GRUT.

**Status: P-24 COMPLETE (read-only census). K1 confirmed as the unique record lock satisfying the census criteria and distinct from the standard classes compared so far (conditional);
K7 flagged unclassified.**

## Addendum 01 — K7 classified; K1 occupancy checked

**K7** (searched the whole repository, including `archive/`, `release/`, `books/`, `uploads/`):

| Quantity | Where it exists | Label |
|---|---|---|
| `R = √(1+α)` | only as founding hypothesis H2 (`founding_h2_R_zeta_bridge`, `to-derive`, FRONTIER-RESERVED): "can R be formulated as a SPECTRAL INVARIANT of a vacuum operator…" — no measurement map | **RESERVED-HYPOTHESIS** |
| `S = 12π/α²` | appears **only as a name** in the `rung9a` ledger note and `BUILD_BANKING_PROMPT.md`; no definition, no observable, no computation anywhere in the repository | **UNIDENTIFIED** |
| `Ω_Λ(α)` | appears **only as a name** in the same two places; no formula. The register's `lambda_undetermined` node states GRUT does **not** fix Λ, so an α-determined Ω_Λ conflicts with the record's own open-field node | **UNIDENTIFIED** (and in tension with `lambda_undetermined`) |
| `S = 12π/(R²−1)²` (α eliminated) | relation between a reserved-hypothesis quantity and an undefined one | **INTERNAL-RELATION**; not a prediction |

**K7 adds no empirical lock.** K1 remains the unique census lock.

**K1** (`probes/K1_OCCUPANCY_01.md`):
- **The z = 0 amplitude form is POINT-INTERSECTION**: Linder's Only Run Gravity (`α_B = 0`) lies exactly
  on `Σ₀ = μ₀/2` at its normalization epoch for every amplitude.
- **The time-dependent form `μη ≡ 1` is DISTINCT-LINE within the declared comparison set**: Only Run
  drifts off by `ε(a)/2`, and leading-order `α_T = 0` Horndeski cannot stay on the line without
  reducing to GR.

**Updated K1 wording:** *unique record lock satisfying the census criteria; its `z = 0` amplitude form is
shared by Only Run Gravity; its time-dependent form (`μη ≡ 1`) has no occupant in the declared
comparison set (beyond-Horndeski/DHOST/nonlocal not evaluated).*

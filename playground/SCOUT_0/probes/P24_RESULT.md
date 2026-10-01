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

> **Census-confirmed (search grade): K1 `Σ₀ = μ₀/2` is the only lock on the record that is
> simultaneously law-level, parameter-eliminated, tied to an identified empirical observable, and off
> a standard class.**
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

**Status: P-24 COMPLETE (read-only census). K1 confirmed as the unique distinctive conditional lock;
K7 flagged unclassified.**

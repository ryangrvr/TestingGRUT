# SCOUT_0 W2 P-23 RESULT — is any law-level invariant left after every supplied input varies?

**Charter:** `PROBE_CHARTERS.md` P-23 (decoherence-rate scaling), widened by the auditor:

> After all supplied GRUT inputs are allowed to vary, is there ANY law-level observable invariant that
> stays fixed and differs from standard open-system physics?

**Governing principle (P-17/P-18):** do not seek to tell deterministic from stochastic ontology using
the same reduced law. Seek discrimination **between allowed reduced laws**.
**Script:** `p23_invariant_checks.py` (log `p23_invariant_checks.log`). **Sources:**
`BASELINE_MAP.md` Tables 1–4; `NO_GO_LEDGER.md`; `PREDICTION_UNIQUENESS_MAP_01.md`;
`GRUT_PREDICTION_GATE_GAMMA_T.md`; `X_FLOOR_MAP.md`; `calc/sigma0_anomaly_screen.py`.
NEW HYPOTHESIS CLASS discipline. Nothing here alters a terminal. The gravity/Π₀ routes are fenced
("NOT REOPENED"): they are read, never computed.

## 0. Verdict

> **In the layer the playground can touch (Level-0 core + the S2/S5/S6 open-system extensions), nothing
> survives.** Every law-level feature that stays fixed while the supplied inputs vary is standard
> open-system structure. Every feature that would differ moves with a supplied input. This is banked as
> a **candidate NON-PREDICTIVITY statement, enumeration grade** (§3).
>
> **One object of the right type exists on the whole record:** the interior-family lock
> **`Σ₀ = μ₀/2`** (equivalently `μη = 1`), with a time-constant shape. Both `x` and `α` are eliminated.
> It is off the conformal scalar–tensor/f(R) line `Σ₀ = 0`. But it rests on a **CHOSEN** sector
> assignment (`p_tt_ansatz`), and GRUT-as-written sets `x = 0`, where it is vacuous. It is a
> **conditional falsifier, not a prediction.** It is the template for what a GRUT prediction must be.

## 1. Candidates, varied against every supplied input

"Moves?" = does the feature change when a supplied item varies (temperature, coupling, bath width or
class, geometry, preparation, coarse-graining, correlation scale, localization, state, `K_b`, drift)?
"Standard?" = is it already a property of standard open-system physics or EFT?

| # | Candidate law-level feature | Moves? | Standard? | Gate |
|---|---|---|---|---|
| L1 | Retained response kernel CM (Level-0 G-D) | no, for any symmetric PD `K` (B2) | **yes**: reversible/Debye relaxation, `Σu²e^{−λt}` | R4 fail |
| L2 | Passivity ⇒ positivity (E-1) | no | **yes** (passive linear response) | R4 fail |
| L3 | Cycle affinity ⇒ ¬CM (E-4) | no | **yes**: broken detailed balance gives complex modes (B2: 3-ring eigenvalues `3 ± 1.386i`) | R4 fail |
| L4 | Noise–dissipation ratio (FDT/KMS) | no | **yes** (BORROWED, record's own label) | R4 fail |
| L5 | C-B discriminator / variance relation | relation fixed; value moves with `T₁` | **yes** (P-23a: generic small-noise expansion) | **retired** (P-23a) |
| L6 | Short-time order of `S − D` (`t², t⁵, t⁶`) | **moves** with bath regime (inertia, correlation time) | yes | R3 fail |
| L7 | Decoherence mass exponent `Γ ∝ m^n` (the charter's object) | **moves**: `n = 2.000` (`R ≪ ℓ`) → `1.84/1.56` (`ℓ = 2`) → `1.04/1.02` (`R ≫ ℓ`) (B3) | yes: any linear mass-density coupling (DP, CSL, gravitational-environment decoherence) | R3 + R4 fail |
| L7′ | Amplitude-free ratio `Γ(m₁)/Γ(m₂) = (m₁/m₂)²` (coherent regime) | holds only for `R ≪ ℓ` (supplied) | **yes**, shared by DP/CSL | R4 fail |
| L8 | Energy-basis decoherence wedge | magnitude needs a staked amplitude | qualitatively distinct from DP/CSL | R5 fail: INVISIBLE-BY-SUPPRESSION, 7–47 orders (NO_GO 4) |
| L9 | Tensor friction `Γ_T ∝ ω³`, fixed coefficient, KMS companion | shape fixed | shape is dimensional analysis: a 4D stress-tensor spectral weight `∝ c·ω⁴` | R5 fail, ≳60 orders (canonical gate §9) |
| L10 | No phantom-divide crossing, `w` one-signed (passive single relaxor) | conditional on one mode and on open rung3 | shared with minimally coupled single-field quintessence and ΛCDM (Vikman 2005) | law-level forbidden region, **non-distinctive**; DESI DR3 exposure stands |
| **L11** | **`Σ₀ = μ₀/2`, `μη = 1`, time-constant** (interior family `μ = 1+xα`, `η = 1/(1+xα)`) | **no**: `x` and `α` both eliminated (B1, exact) | **no** for the conformal class (f(R), DGP: `Σ₀ = 0`); whether some other standard theory lies on the line is **not established here** | R1 **fail** (rests on CHOSEN `p_tt_ansatz`; family "inherited-conditional"); R3 pass (the relation is parameter-free); R4 pass vs conformal class; R5 open (testable only if `x ≠ 0`; no floor on `x` — no-pin) |
| L12 | Matsubara comb at spacing exactly `H` (prediction map P2) | — | — | **fenced** (gravity sector, not reopened); unrun upstream |

**Explicit parameter eliminations derived:**

- L11: `Σ − 1 = (μ − 1)/2`, independent of `x` and `α` (B1).
- L7′: amplitude cancels in `Γ(m₁)/Γ(m₂)`, but the exponent is fixed only by supplied geometry.
- L5: bath amplitude and width cancel against measured `v` (P-23a).

**Only L11 survives the variation and differs from a standard class.**

## 2. Reading L11 honestly

- **Kind:** a cross-observable, parameter-eliminated, sign-and-ratio relation. It is exactly the kind of
  object the auditor's list asks for (cross-observable ratio; forbidden region off the line).
- **Already on record:** `calc/sigma0_anomaly_screen.py` calls it the family's "EXACT LOCK,
  inherited-conditional". **REDISCOVERED-RECORD.** The scout adds the classification: it is the
  **only** surviving law-level invariant across the scanned layers.
- **Why it is not a prediction:** (i) it inherits the CHOSEN projector assignment (R1); (ii) at GRUT's
  constitutive choice `x = 0` it reduces to ΛCDM and tests nothing; (iii) `x` has a ceiling but no
  floor, so a null is uninformative.
- **What it can do:** act as a **conditional falsifier** of the interior family. A resolved
  `(μ₀, Σ₀)` off the line, or a late-switching (`Ω_DE`-tracking) shape instead of a time-constant
  one, excludes the family. A point on the line would not confirm GRUT unless the line is shown to be
  unoccupied by other theories (an owed literature check, not done here).
- **What would make it a prediction:** deriving the sector assignment instead of choosing it, plus a
  floor on `x`. That is the canonical Π₀ route (X_FLOOR D2b), which is fenced.

## 3. Candidate NON-PREDICTIVITY statement (enumeration grade)

> **Within the earned Level-0 core and its declared S2/S5/S6 open-system extensions, the set of reduced
> laws reachable by varying the supplied inputs is restricted only by:**
> - **(a)** reversible linear-response structure (symmetric PD coupling ⇒ CM, positivity);
> - **(b)** FDT/KMS;
> - **(c)** broken-detailed-balance signatures (affinity ⇒ ¬CM);
> - **(d)** diffusion-limit universality (P-17/P-18: only `Q` survives).
>
> **All four are shared with standard open-system physics. So this layer does not restrict the reduced
> law beyond structures already shared by standard effective theories.**

- **Basis:** enumeration of `BASELINE_MAP.md` Table 1. Each E-entry is one of:
  - a standard structural property (E-1, E-2, E-4, E-7, E-9, E-11);
  - a certified non-implication or terminal about selection (E-5, E-6, E-8, E-10, E-18, E-20);
  - a class result conditional on supplied inputs (E-12, E-13, E-14, E-15, E-16, E-17, E-19, E-21,
    E-22).

  None fixes a law-level quantity outside (a)–(d). Table 2's conditional results inherit their
  supplied inputs by definition.
- **Grade:** enumeration, not proof. It holds if the Table 1 enumeration is complete and correctly
  classified; one mis-classified entry would break it. **Strongest objection:** "standard" is judged
  per entry, not by one formal criterion. A formal version would need a definition of "standard
  effective theory class" (e.g. Onsager–Casimir linear response + FDT + Markovian/GLE reductions) and
  a proof that the earned predicates generate nothing outside it.
- **Scope fence:** the cosmology/gravity sector is **not** covered. There, L11 is the one live
  law-level object, plus the fenced Π₀ and Matsubara-comb routes.

## 4. What this changes

- **Target sharpened (auditor's reduction, now evidenced):** GRUT cannot gain empirical
  distinctiveness by assigning a different hidden ontology to the same effective law (P-17/P-18). In
  the open-system layer it **also** does not restrict the law beyond standard structure (§3). Any
  distinctiveness must come from a sector that fixes a law-level relation. On the current record, the
  only such relation is L11, which is conditional on a choice.
- **For the program (owner's call, not opened here):** the decisive in-record question becomes
  whether the `p_tt` sector assignment, and with it the `Σ₀ = μ₀/2` lock, can be **derived**. That is
  the X_FLOOR / Π₀ territory. The scout does not enter it.

**Status: P-23 COMPLETE (one pass). Open-system layer: nothing survives → candidate NON-PREDICTIVITY
statement (enumeration grade). Whole record: exactly one law-level invariant of the right type
(`Σ₀ = μ₀/2`), conditional on a CHOSEN projector → conditional falsifier, not a prediction.
Decoherence exponent (charter object) moves with supplied geometry and is shared with DP/CSL.**

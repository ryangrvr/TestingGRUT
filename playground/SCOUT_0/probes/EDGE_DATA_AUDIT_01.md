# SCOUT_0 EDGE-DATA AUDIT 01 (read-only; requested after ZOOM_OUT_03)

**Question:** for every earned result E-1…E-22 (`BASELINE_MAP.md` Table 1), does it **independently
FIX** a local spectral edge datum relevant to an effective law, or does it only state a consequence
once that datum is supplied?

**Edge data audited:**
- edge location / spectral bottom;
- edge exponent (local density-of-states exponent);
- edge velocity at an occupation edge;
- visible spectral sign;
- number and location of soft points.

**Equivalent invariant found during the audit:** the **low spectral moments** of `μ_r` (e.g.
`K₁₁ = ∫λ dμ_r`), which control **short-time** coefficients. Edge data control the IR.

**Witness script:** `eda_witnesses.py` (log `eda_witnesses.log`). Read-only otherwise. Candidate
structure only; no terminal touched.

## 0. Verdict

> **No E-entry is FIXES-EDGE-DATA.** Every law-level predicate the earned layer certifies either
> (a) **reads out** edge or moment data of a local spectral measure, or (b) is independent of
> spectral data (selection terminals, non-implications, architecture).
>
> Every edge or moment **value** traces upstream to a supplied primitive: the coupling matrix/graph
> `K` (S-1), the pin, the sector (S-7), statistics/exclusion (S-7, P-02), the temperature/noise
> profile (S-8/L0-1e), the drift (S-9), the preparation (S-11), the readout site (S-10), or the
> gravitational inputs (S-6, S-13).
>
> Explicit two-model witnesses exist (§3): same earned predicates, different supplied primitive,
> different edge datum, different effective law. **No earned edge selector was found**, so P-02b is
> opened next as directed.

## 1. Edge-dependency matrix

Classes: **FIX** = FIXES-EDGE-DATA · **CON** = CONSTRAINS-BUT-DOES-NOT-FIX · **READ** =
READS-OUT-EDGE-DATA · **IND** = INDEPENDENT · **SUP** = SUPPLIED-UPSTREAM · **NF** = NOT-FORMULABLE.

| E | Earned result (abridged) | Class | Edge/spectral datum involved | Traced upstream to |
|---|---|---|---|---|
| E-1 | gap → P_memory; passivity → P_positivity | **READ** | spectral bottom; visible sign (P-01: both are `μ_r` statements) | pin, `K` (S-1). Witness W-A |
| E-2 | locality → P_geometry | **CON** | locality makes `K` banded/finite-range (compact `supp μ_r`), but bottom and exponent stay free: chain vs half-plane (W-A) and chain vs expander (P-01) are all local | graph/dimension in `K` (S-1) |
| E-3 | linearity → exact reduction (definitional) | **CON** | makes the spectral representation `k = ∫e^{−λτ}dμ_r` exist at all; fixes no value | — (definitional) |
| E-4 | cycle affinity → ¬P^resp_CM | **READ** | visible sign / complex modes | asymmetric `K` (supplied affinity) |
| E-5 | seven certified non-implications | **IND** | several are themselves **non-fixing** statements (e.g. spectral stability ⇏ accretivity ⇏ monotone response); they support the verdict | — |
| E-6 | no single deep organizing pair | **IND** | — | — |
| E-7 | TRIVIAL/IDENTITY observability at end-site readout | **IND** | makes the retained object site-relative (`μ_r` at the declared site) | readout site (S-10) |
| E-8 | P-1/X2 intrinsic partition selection (one regime, one of three criteria) | **CON** | the one earned **member selection**: it picks *which* site is special (the impurity), so it can fix *which* `μ_r` is read. The impurity, its strength and `K` are in the supplied model, so the edge values still follow from supplied data | supplied substrate containing the impurity (S-1) |
| E-9 | earned spine; P-2 influence cone CONE-CONFIRMED (Gaussian class) | **READ** | cone speed = maximal group velocity of the dispersion (a velocity datum) | `K`, supplied ħ/quantum class (S-4, C-6) |
| E-10 | D-1: medium/relational vocabularies demoted | **IND** | — | — |
| E-11 | exact FRW kernel transport; clock objection retired | **SUP** | transports a kernel whose spectral data are inputs | S-13 (FRW state, KMS/T, α, …) |
| E-12 | GR2-L6: finite-window counting asymptote `S_∞ ≈ 6.2142 = 2d + 0.214` | **READ** | density-of-states / counting exponent (`2d`) plus a lattice-dispersion correction | dimension, lattice, window (declared counting class) |
| E-13 | SF-1 FORMATION-OF-LAW-CLASS | **READ** | edge velocity at the occupation edge; soft-point count (P-02 theorem) | sector (S-7) + exclusion (S-7; P-02) |
| E-14 | S5-1 MARKOV-LIMIT-OTHER-CLASS; native `t^{−3/2}` branch cut | **READ** | `t^{−3/2}` = chain edge exponent `γ = 1/2` (ZOOM_OUT_03 check); `κ = 1/(2√K₁₁)` = first moment | `K_b`, pin (W-C: κ moves 0.3297 → 0.3450 / 0.3162 as pin 0.3 → 0.1 / 0.5) |
| E-15 | S2-1 FULL-DISCRIMINATOR-CONFIRMED | **READ** (moments; nonlinear class) | `Δc₂ = −24βT₁a` reads `T₁`; `Δc₃` reads `K₁₁ = ∫λ dμ_r` (W-B); at O(β) the discriminator is a functional of linear-response `g`, `v` (P-17) | drift β (S-9), noise profile (S-8), `K` (S-1) |
| E-16 | S2-HB UNRESTRICTED-REALIZATION-ONLY | **IND** | realization, not law (P-17/P-18) | — |
| E-17 | S3-0 FORMULABLE-ONLY-WITH-CHANGE | **NF** | crossed cell not in the declared theory | — |
| E-18 | REV-0 NOT SUPPORTED | **IND** | — | — |
| E-19 | S6 net arrow; `D, Ḋ ~ t^{−6}` tails | **READ** | tail exponent from the parent's band edges | parent, preparation (S-11), J bookkeeping |
| E-20 | SYN-0 architecture | **IND** | meta | — |
| E-21 | sector-selection S-1 PROVISIONAL (restated: counting does not select the exponent) | **SUP** | "amplitude and state sector-supplied" | sector (S-7) |
| E-22 | CA-1 earned carrier elimination (within supplied quantum class) | **IND** / **SUP** | carrier, not edge | supplied quantum class (A-6) |

**Tally:** FIX 0 · CON 3 (E-2, E-3, E-8) · READ 8 (E-1, E-4, E-9, E-12, E-13, E-14, E-15, E-19) ·
IND 8 (E-5, E-6, E-7, E-10, E-16, E-18, E-20, E-22) · SUP 2 (E-11, E-21) · NF 1 (E-17).

## 2. Dependency test (upstream trace)

Every READ/CON row ends at a supplied primitive, never at an earned statement that fixes a value:

| Primitive | Edge/moment datum it controls | Earned theorems stay true when it varies? |
|---|---|---|
| pin | spectral bottom | yes (E-1 is an implication; P-01) |
| graph / coupling `K` | bottom, exponent, velocity, moments | yes (W-A, P-01) |
| sector | occupation-edge location → edge velocity | yes (SF-1 itself) |
| statistics / exclusion | map ν ↦ `k_b^∞` | yes (P-02) |
| noise profile `T_i` | the amplitude read by Δc₂ | yes (S2-1 identity holds for every profile) |
| drift β | curvature read by C-B coefficients | yes |
| preparation | which edge contributes (S6 tails) | yes |
| readout site | which `μ_r` | yes (E-7; E-8 selects only within a supplied impurity model) |
| complex/quantum lift | — (P-08/P-09: observations forget it) | yes |
| gravitational inputs | kernel spectral data in transport | yes (E-11) |

## 3. Two-model witnesses (same earned predicates → different edge → different law)

- **W-A (graph):** pinned chain (n = 400) vs pinned 2D half-lattice (46×46), pin = 0.3, both
  nearest-neighbour.
  - Both satisfy every E-1/E-2/E-3 predicate: local, passive, gapped (`λ_min = 0.3000 / 0.3002`),
    CM kernel.
  - The pin-stripped kernel `k₀ = e^{pin·τ}k` has local slope **−1.471 → −1.497** (chain:
    `γ = 1/2`, `τ^{−3/2}`) vs **−1.18 → −1.15** (half-plane boundary site: 2D-type `τ^{−1}` with a
    slowly drifting logarithmic point-defect correction).
  - Different edge exponent, different memory law (prefactor universality class), same earned
    predicates.
- **W-B (moments):** `∫λ dμ_r = K₁₁` exactly (2.300000 chain, 4.300000 half-plane). The record's
  C-B `Δc₃ = 24βT₁a(44βa² + 5K₁₁)` therefore reads a supplied first moment.
- **W-C (value):** S5-1's `κ = 1/(2√K₁₁)` moves 0.3297 → 0.3450 / 0.3162 as the supplied pin goes
  0.3 → 0.1 / 0.5. The class is fixed (Markov limit); the value is not.
- **From earlier probes:**
  - SF-1 sectors D vs E: same parent, different occupation edge, z = 1 vs 2.
  - P-02 statistics: same sector, exclusion vs free bosons, (1,2) vs (2,1).
  - P-01: chain vs expander at pin 0, algebraic vs exponential memory.

## 4. Theorem candidate (scoped as instructed)

> **Within the currently certified quadratic/free GRUT constructions** (the C1-a/L0-1 chain family
> and its declared linear members, the SF-1 free parent, S5-1's linear core, the S6 parent, and the
> linear-response layer of S2), **all earned law-class predicates factor through local spectral data
> of a retained observable** — edge data (location, exponent, velocity, visible sign, soft points)
> for IR predicates, low moments for short-time predicates. **The values of those data are functions
> of supplied upstream structure.** Therefore the earned layer classifies the consequences of
> spectral data but does not independently select the effective law class.

- **Grade:** audit + witnesses (each READ row's dependence is exhibited, and four rows have explicit
  witnesses). Not a proof over all conceivable earned statements.
- **Not claimed:** anything about interacting/non-quadratic parents; the global iff "GRUT restricts the
  law iff it fixes edge data" stays a **conjecture** until P-02b (or another interacting parent)
  tests whether spectral data remain the complete quotient.

**Strongest obstruction:**

1. **E-8** is a genuine earned member selection. In a richer model where the selected site's `μ_r`
   were fixed by the selection criterion itself (not by supplied impurity data), it could become a
   CON → FIX upgrade. The audit found no such case at the recorded scope.
2. **E-15/S2 is nonlinear.** Its edge/moment reading holds at O(β) through linear response. Its
   higher orders (`β²`, `κ₄` — P-17 C4) involve data beyond `μ_r` (drift nonlinearity, noise cumulants),
   which are supplied but are not spectral data. The quotient there is (`μ_r`, drift, noise law), not
   `μ_r` alone.
3. "Factor through" is shown row by row, not from one formal definition of "law-class predicate".

## 5. Relation to P-23

- **STRUCTURAL EXPLANATION OF THE ENUMERATION: YES.** P-23 §3 listed four kinds of fixed features:
  - (a) CM/positivity of reversible response = visible sign + real spectrum of `μ_r`;
  - (b) FDT/KMS = supplied noise rule;
  - (c) affinity ⇒ ¬CM = complex/visible-sign data;
  - (d) diffusion universality = only the noise intensity `Q` survives.

  Items (a) and (c) are exactly the "READS-OUT" rows. (b) and (d) are supplied noise-side data. The
  enumeration now has a mechanism: **the earned layer certifies functions of spectral data and never
  fixes the data.**
- **PROOF OF GLOBAL NON-PREDICTIVITY: NO.** That would need every law-level invariant outside the
  audited class to be excluded. The audit cannot do that: interacting parents, gravity-sector kernels
  and the fenced Π₀ route are untested.
- **Grade change for P-23 §3:** from "enumeration grade" to **"enumeration + mechanism (audit
  grade)"** within the quadratic/free scope. The statement is not upgraded beyond that scope.

**Decision:** no earned edge selector found → **P-02b opened next** (interacting bosons; the free-boson
P-02 result kept as the unmodified control).

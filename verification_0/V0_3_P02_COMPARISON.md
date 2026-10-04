# V0-3 — COMPARISON PHASE (P-02 / HCB-SF1; the original was unsealed only after the reproduction was committed at `4ce8686`)

**Original:**
- `playground/SCOUT_0/probes/P02_RESULT.md` (blob `f766a727…`);
- `P02_SF1_CONTROL_MAP.md`;
- `p02_boson_sectors.py` (definition and print lines only, for comparison).

**Reproduction:** `V0_3_P02_REPRODUCTION.md` and `code/v0_3_*`. These were produced by a context-isolated sub-agent from
the spec alone. They are **not** edited below.

## 1. Scientific comparison

| item | original | independent reproduction | classification |
|---|---|---|---|
| Free-boson ground state | unique condensate (b₀†)^N\|0⟩/√N!, E₀ = Nε(0) | identical (unique; E₀ = −2N) | **exact agreement** |
| ρ_q action / support | √N\|(N−1)₀, 1_q⟩, single line 4sin²(q/2), weight N/N = 1 | identical for q ≠ 0 | **agreement** |
| **"S independent of N"** | at every L (unqualified) | holds for **q ≠ 0 only**: S(q = 0) = N·δ(ω) | **precision (CR-a)**. q = 0 is outside P's q ∈ (0, π], so the class is unaffected |
| Zero-mode treatment | k = 0 is the mechanism, not an artefact; no regulator needed | the same; additionally requires a **unique and quadratic** minimum | **agreement plus scope (CR-b)**. The original's hostile set did not test degenerate or flat minima. With a quartic minimum, bosons give (4, 1); with degenerate minima the class depends on the declared reference state |
| Collapse class | (2, 1) in every family, C-6 included, under P and P′ | identical for the cosine band | **exact agreement** |
| A1 (L ∈ {6, 10}, N ∈ {1, 3, 5}) | unique GS; gap 1.000 / 0.382; single line | unique; gap 1.000000 / 0.381966; single line | **exact agreement** |
| HCB ↔ SF-1 (JW) | twist (−1)^{N−1}; periodic for odd N; ρ JW-invariant; exact for all odd N at L = 6, 10, 12; 1D-scoped | identical, 14 / 14 odd-N sectors; even N fails (L = 12, N = 6: E₀ −7.7274 vs −7.4641); **also needs nearest-neighbour hopping and a periodic ring**; d > 1 fails (4×4 torus) | **agreement plus scope (CR-c)**: the NN-hopping / periodic-ring conditions are added to the original's 1D / odd-N scope |
| "Exclusion, not antisymmetry, is the carrier" | stated, 1D-scoped | confirmed within those conditions | **agreement** |
| SF-1 classes and v_b (A3) | D, D-¼, D-¾ = (1, 2), v_b = 2, √2, √2; E, E-3, Ē = (2, 1), v_b = 0; z_P = z_{P′} | identical | **exact agreement** |
| C-6 | k_b ~ L^{−1/2}; the two-scale approach of v_b → 0 competes with q ~ 1/L; path-dependent | identical mechanism: z_P = 2, z_{P′} = 3/2; the path exponent is 2 for a ≤ ½ and 1 + 1/(2a) for a > ½; riding the soft point gives 3 | **agreement** (the reproduction quantifies the path dependence) |
| A4 (Fermi pattern in bosons) | E − E_gs = 372.8 at L = 1026; min ω = −1.99, −0.76, −0.20 | sector identified as D, N = 513: 372.8271; −1.993876, −0.763484, −0.195623 | **exact agreement** |
| **Exponent rule** | "z = 1 iff ε′(k_b^∞) ≠ 0; **z = 2 iff k_b^∞ is a band extremum**" (§0, §4). §5 adds "(quadratic extremum)" | **false in general:** z = r_e, the order of the first non-vanishing derivative at the edge. Counterexamples: quartic minimum gives z = 4; stationary inflection gives z = 3 | **CORRECTION (CR-d)**: "extremum ⇒ z = 2" holds only for quadratic extrema |
| **"Law class = function of the single invariant k_b^∞ (equivalently v_b)"** | stated as a theorem for **free, number-conserving, single-band, smooth ε on a 1D ring** | **refuted within that declared class:**<br>• v_b does not fix the soft count: v_b = 0 occurs with I-q of 1 or 2, and v_b ≠ 0 with 2, 3 or 10;<br>• multiple pockets (ε = −2cos3k − ½cos k) give I-q = 10;<br>• for asymmetric bands, k_b = πν fails and P sees edge orientation;<br>• degenerate minima make the class reference-state dependent. | **CORRECTION (CR-e), grade-bearing.** The theorem holds only under extra assumptions: an even band, ε′ ≠ 0 on (0, π), quadratic extrema, a single contiguous Fermi sea, fermions at fixed ν (the cosine band qualifies). Even then, only **whether** v_b = 0 matters for z. **Corrected invariant:** the set of edge data (positions, orientation, local order r_e), with z = max r_e over the edges P probes, and Q_soft = the set of pairwise edge differences modulo the quotient. **Exponent and soft count are separate invariants** |
| Two-band case | outside the declared single-band class | two-leg ladder gives I-q = 3 with two velocities | outside the original's class: **scope note only** |
| "Statistics act only through ν ↦ k_b" | free bosons: k_b ≡ argmin ε (pinned) | for bosons, ν plays no role at all; the condensate "occupation edge" is shared notation, not the same object as a Fermi edge (they coincide only as N = 1 / ν → 0) | **precision (CR-f)**: compatible with the original if "the map ν ↦ k_b is constant" is read literally; the conceptual distinction should be recorded |
| Fourier convention | ρ_q = Σ_j e^{−iqj}n_j alongside Σ b†_{k+q}b_k | the two imply opposite ε(k) orientation conventions; irrelevant for even bands, relevant for asymmetric ones | **precision (CR-g)** |
| Interaction conjecture (P-02b) | spawned; not part of P-02 | not used (firewall) | out of scope |

## 2. Grade

**V0-3-C — P-02 / HCB-SF1 REPRODUCED WITH CORRECTION — VER-I1 (orchestrator-exposed).**

**Core results reproduced exactly:**
- free-boson collapse to (2, 1) in every family under P and P′;
- the unique condensate;
- the single excitation line;
- A1 – A4 numbers;
- HCB ↔ SF-1 exact equality for all 14 odd-N sectors at L = 6, 10, 12, via Jordan–Wigner with twist (−1)^{N−1} and a
  JW-invariant density;
- exclusion, not exchange antisymmetry, as the 1D carrier;
- the cosine-band SF-1 classes and v_b values;
- the C-6 two-scale mechanism.

**Corrections:**
- **CR-e (grade-bearing; the theorem's generality is refuted within its declared class):** "the law class is a function of
  the single sector invariant k_b^∞ / v_b" is **false** for general free single-band 1D parents. **Banked corrected
  statement:**
  - **Exponent:** z = max r_e over the edges probed by P, where r_e is the order of the first non-zero derivative of ε.
  - **Soft count:** I-q = the number of inequivalent pairwise edge differences, modulo q ~ −q ~ q + 2π.
  - These are separate invariants of the **edge-data set** (positions, orientation, local jets).
  - The single-scalar form holds only for even, generic bands with a contiguous Fermi sea, such as the cosine band. There
    only the vanishing or non-vanishing of v_b matters, not its value.
- **CR-d:** "z = 2 iff extremum" becomes "z = 2 iff the edge is a **quadratic** extremum"; in general z = r_e.
- **CR-b:** boson collapse to (2, 1) needs a **unique quadratic** minimum. In general it is (r\*, 1); for degenerate
  minima the class is undefined until a reference state is declared.
- **CR-c:** the HCB ↔ SF-1 match needs 1D, **nearest-neighbour hopping**, a **periodic ring** and odd N.
- **CR-a, CR-f, CR-g:** precisions (S at q = 0; the meaning of "occupation edge" for condensates; the Fourier
  convention).

**Why C and not F.** Every core computational and mechanistic result (collapse, HCB match, carrier, SF-1 classes) is
reproduced exactly. What fails is the **generality** of the stated sector-invariant theorem. The rubric places that under
C ("theorem generality, invariant definition, exponent statement, soft-count statement … requires correction"). **Owner
option:** if the single-v_b theorem is regarded as a core claim of P-02, the TARGET-3 sub-result would merit F,
**refuted as stated, with corrected replacement**, while the rest stays C.

**Independence:** VER-I1, orchestrator-exposed. The orchestrator saw the original's short proof lines during extraction.
The reproducer was context-isolated, spec-only, with a file-access report. **Not VER-I2.**

**Frozen `scout-0` is not modified.** The corrections live in VER0 and in later syntheses citing P-02. The original's §6
proposal ("the law class factors through v_b^∞") should be read with CR-e.

## 3. Owner ruling (additive; review of `9d673bb`)

**V0-3 ACCEPTED.**
- **Overall grade:** **V0-3-C — P-02 / HCB-SF1 REPRODUCED WITH CORRECTION — VER-I1 (orchestrator-exposed).**
- **Explicit subgrade:** **V0-3-T3-F — SINGLE-v_b SECTOR-INVARIANT THEOREM REFUTED AS STATED.**
- **Criterion-2 item 2:** COMPLETE AT VER-I1 WITH CORRECTION.

**Survives:**
- the free-boson collapse;
- HCB ↔ SF-1, at its corrected scope;
- exclusion as the 1D carrier;
- the registered cosine-band sector results;
- A1 – A4;
- C-6.

**False within its declared class:** the theorem that the general free 1D single-band law class factors through one
scalar k_b^∞ / v_b.

**Banked replacement.** The relevant quotient is an **edge-data set**: positions and orientations of the relevant
occupation edges, plus their local dispersion orders / jets. The dynamical exponent and the soft-point structure are
separate invariants. At the verified scope:
- local support exponents depend on r_e at the relevant edges;
- soft momenta depend on the inequivalent allowed edge differences;
- multiple pockets cannot be reduced to one v_b;
- the scalar v_b = 0 / ≠ 0 classification is only a special compression for the generic single-pocket cosine-like case.

**CR-a … CR-g are accepted.** Frozen `scout-0` is not modified.

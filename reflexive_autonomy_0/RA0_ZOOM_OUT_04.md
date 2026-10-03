# RA0 ZOOM-OUT 04 (after RA0 REPAIR 02 + G4) — HARD STOP for owner review

**Branch:** `grut-reflexive-autonomy-0`. New theory construction. Not frozen. No consciousness, quantum, collapse or
Born-rule work.

## 1. Does the canonical slow projector asymptotically close under multiplication?

**Yes, in the nearly-uncoupled reversible class (PROP G4-T, assembled from Weyl and Davis–Kahan):**
- Δ_alg ≤ s·(4/√p_min + C_V), with s ≤ η/(g − η) → 0.
- Logged values:
  - F1: Δ_sup 0.31 → 0.072 (∝ 1/M);
  - F2: ∝ M^{−2} at the coarse level and ∝ M^{−1} at the fine level;
  - modulated ladder β = 3: 7.6e-4 → 3.8e-6.

**No in SSEP, independent particles and diffusive ladders:** forced slow spaces keep Δ_sup ≈ 0.7 – 1.3, flat in size.

## 2. Is there a theorem linking approximate closure to a nearby exact partition algebra?

**Partly.**
- **Rank 2, PROVED HERE (PROP G4-P), dimension-free:** ‖E_V − E_{A_f}‖ ≤ (κ − 1 − s²)^{1/2} ≤ Δ_alg, where
  A_f = {f > s/2}.
- **Nearly-uncoupled class:** Davis–Kahan gives a nearby partition algebra.
- **General rank ≥ 3: CONJECTURE G4-C** (gap ≤ c_k·Δ_alg). It held in 40/40 logged runs (max ratio 0.93), but is not
  proved.
- **Literature:** no general theorem was located. The closest is the Peng–Sun–Zanetti structure theorem, but its
  hypothesis is λ_{k+1}/ρ(k), not λ_{k+1}/λ_k. AMNM and Kadison–Kastler concern maps or pairs of algebras, not almost-closed
  subspaces.

## 3. Can the partition be recovered blindly from P_N?

**Yes: BLIND PARTITION RECOVERY.**
- **Procedure:** the spectral data of P_N, π (derived), and the certified cut rank give attracting tensor-power fixed
  points, which are canonical idempotents; argmax over them gives the partition.
- **Results:**
  - The number of blocks is an output, and it equals the cut rank in every positive run.
  - Exact match to the hidden labels, which were checked only afterwards:
    - F1 for M ≥ 20;
    - F2 super-basins at all M, and F2 basins for M ≥ 10;
    - ladder lanes at all L.
  - F1 at M = 5 and 10 fails. That regime is pre-asymptotic (η > g; no rank-3 cut has formed yet).
- **Proof status:** proved for rank 2; numerical for rank ≥ 3.

## 4. Does the procedure reject SSEP / independent-particle false positives?

**Yes.** SSEP, independent particles and ladders with β ≤ 2 are refused, because none has a certified diverging cut
(PROP G3-A, PROP G3-B, bounded ladder ratios). Forced stress runs confirm that there is no algebra in SSEP, independent
particles or diffusive ladders.

**Two honest findings:**
- **Ladder β = 2:** the cut is bounded (≈ 19.4), so it is refused, yet the forced slow space **does** algebraize (lanes,
  Δ → 7.7e-4). Algebraization and a diverging cut are distinct properties; the cut rule is conservative.
- **Finite-size ratio trends cannot certify divergence:** SSEP's ratio rises 1.62 → 1.78 toward a bounded limit, while
  F1's ratio dips before diverging. Certification must be analytic.

## 5. Does a nested spectral filtration induce a nested partition filtration?

**Yes, numerically (F2).** The blind rank-2 partition is coarser than the blind rank-4 partition, with **zero violating
π-mass at every M**, in the order of the spectral filtration. This matches the hidden super-basins and basins.

**Terminal: CANONICAL DYNAMICAL PARTITION FILTRATION.** No level is selected.

## 6. Did any threshold, K, labels or coordinates enter?

**No:**
- labels and coordinates are used only by the post-hoc audit, which is separated in the code;
- K is an output;
- there is no modelling threshold (computational equality tolerances only);
- there is no crispness objective.

**Priced:**
- (a) the rank comes from an **analytic certificate of divergence, which is a proof about the supplied family**;
- (b) the L²(π) norm (reversibility, irreducibility);
- (c) the family itself.

## 7. Does the result derive A_partition in the metastable class, or only a timescale sector?

**More than a timescale sector: dynamics → observable algebra → partition**, conditionally, in class G4-T.

But the partition is a **coarse-graining of states** (macrostates), modulo vanishing π-mass. It is a derived resolution /
coarse-graining **given the family**.

It is **not** a subsystem decomposition (tensor factor or net), which is what Σ / A_partition name in frozen GRUT. A
metastable partition does not split the system into interacting parts, so **A_partition is not derived.**

## 8. Is any frozen residual input eliminated?

**No.**
- The coarse-graining is no longer chosen once the family is fixed. But the family (and hence its scale separation) is
  supplied, so A_resolution is **relocated into the family**, as in G2 and G3.
- Σ / A_partition are untouched.
- Orientation: all families are reversible, so nothing is selected.
- No frozen witness is distinguished. **Not TRUE COMPRESSION.**

## 9. Is a consciousness construction finally mathematically eligible?

**Eligible to be *posed*, not justified.**

**What G4 supplies.** For the class G4-T, G4 supplies what G3 lacked: an emergent observable algebra and partition
(filtration) encoded in D, with no supplied boundary, ε or K. The question "can an internally realized process be defined
relative to this emergent algebra, without an externally supplied subject or system boundary?" therefore now has a
mathematical object to refer to.

**What G4 does not supply:**
- awareness, experience, subjectivity, readout, self-model, memory or observer selection;
- a subsystem / tensor-factor structure: the emergent object is a macrostate partition, so any "internal process" would
  have to be defined on macrostate dynamics, not on a part of the system.

**Recommendation, if the owner opens anything:**
- restrict it to class G4-T;
- preregister that the family remains supplied;
- first close CONJECTURE G4-C (rank ≥ 3), or work only at rank 2, where the derivation is fully proved.

## Scorecard

| gate | terminal |
|---|---|
| G1 | FINITE REFLEXIVE AUTONOMY NO-GO |
| G2 | CONDITIONAL MACROVARIABLE DERIVATION FROM A SUPPLIED FAMILY |
| G3 (repaired RA2) | CONDITIONAL DYNAMICAL HIERARCHY DERIVATION |
| G4 | **CONDITIONAL DYNAMICAL PARTITION DERIVATION** (class G4-T); BLIND PARTITION RECOVERY (proved rank 2, numerical rank ≥ 3); CANONICAL DYNAMICAL PARTITION FILTRATION (numerical); false positives refused; algebraization ≠ diverging cut (ladder β = 2) |
| frozen witnesses distinguished | 0 |
| residual inputs eliminated | 0 |

**HARD STOP.** Returning for owner review.

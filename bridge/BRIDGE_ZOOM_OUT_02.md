# BRIDGE ZOOM-OUT 02 (after BRIDGE REPAIR 01 + B3) — HARD STOP for owner review

**Branch:** `grut-bridge-1`. **Canonical GRUT:** read-only at `master-w25bu9 @ b935099`, untouched.
**Evidence:** `B3_ACCESS_RESULT.md`, `B3_ACCESS_COMPONENT_LEDGER.md`, `b3/b3_access.log`. B1 scope repairs are in
`BRIDGE_CORRECTION_LEDGER.md`.

## 1. Did GRUT eliminate any component of A_res?

**No component was eliminated (TRUE COMPRESSION = 0).**

**One component moved downstream conditionally.** The closure / effect algebra is a definite function of (D, seed) once a
closure rule is fixed:
- Classical: Krylov and spectral methods agree to 1e-14.
- Quantum (auxiliary): reproduces P-5's 16 / 8 / 4 ladder.

That is **CONDITIONAL COMPRESSION**. It is the first non-zero Bridge compression, and it is modest, as anticipated.

## 2. Is access best described as one primitive or several?

**Several.** The tested classes need **four supplied blocks plus one convention** (A_interface is itself a pair, seed + readout [BR2-01]):

    A_res  →  A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time,   A_interface := (A_seed, A_readout)   [BR2-01]
    A_closure = f(D, A_seed; R_closure)   [BR2-02]

- **A_partition is not reducible to A_seed.** The same total readout grouped differently changes the record counts.
- **A_coarse ⊆ A_partition** is shown only for fragment-restriction coarse maps (CONDITIONAL; aggregating maps not tested).

## 3. Is A_closure genuinely derived given A_seed?

**Yes, conditionally, and with a caveat the canonical record had not surfaced.**
- **Given a rule**, the closure is unique, reproducible and representation-covariant.
- **The canonical instruments use more than one rule:**
  - P-5 uses alg{I, H, B};
  - P-6 uses an ad_H-closed algebra of the seed;
  - P-6 L-D2 uses no dynamics at all.
- These agree on generic, parity and abelian cases. They **disagree** when H has structure outside the seed's reach:
  dimension **16 vs 4** for a site decoupled from a hidden sector.

Proposed **B3-CUC-1** (bookkeeping clarification, not applied):
- book the closure as downstream of (D, seed);
- standardize one rule (R_P6 is the natural candidate; owner's call).

No REDUNDANT SUPPLY was found. Ledger A-13's "access sets" are declared site sets, not closure algebras.

## 4. Does A_seed remain independent of D + Σ?

**Yes: A_seed NOT REDUCIBLE TO D + Σ.** The no-go witness has the same D, the same Σ, an identical closure and center,
and different downstream physics:
- **Classical earned class:** 112 / 120 distinct pairs among complete single-site seeds on C1-a. All 276 pairs are
  distinct on an asymmetric weighted path, including the two **end sites** (same Σ role).
- **Auxiliary quantum class:** σ_z¹ vs σ_z³ on a parity-conserving 3-spin chain. The closures are identical (32-dim, center
  2), and the probe coherence differs by 0.52–0.67 (3/3 trials). This re-establishes P-6 L-B.

Minimality does not rescue selection:
- backward elimination lands on different minimal sets ({1, 2, 3, 4} vs {0, 1, 2, 3}; five tied minimal sets);
- informational completeness holds for an open dense set of readouts.

## 5. Did EA-0 provide compression, relocation or only an auxiliary-class construction?

**At earned Level-0: BLOCKED / CATEGORY MISMATCH.** EA-0 is UNFORMULABLE, as carried exactly. No ρ or quantum reach was
imported.

**In the auxiliary class:**
- faithful states collapse access to the identity;
- the generic case is non-selective;
- non-trivial access appears only where D has exact block / invariant structure. That is **RELOCATION INTO D**, not
  access elimination.

Theorem S is sufficient only. The retired "generator-carried iff" is **not** revived. EA-0 provided **no compression**.

## 6. Can GS1 be inverted from geometry to access?

**No: GS1 IS ONE-WAY (A → geometry).** On one fixed hidden grid:
- every tested declaration is admissible: full, two different 3-site boundaries, and two single sites;
- **all 9 singletons and all 84 three-sets** are informationally complete;
- each declaration gives different interface data;
- every candidate inverse selector fails: max-information returns the full set, min-cardinality gives 9 ties, and GS1
  selection power selects full access only.

The access boundary stays supplied.

## 7. Did the reviewed residual shrink?

**The boundary did not shrink. The A_res component was sharpened.**

    C5 → D_dyn ⊕ [Σ ⊗ H_corr]_coupled ⊕ A_res,
    A_res = A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time   (A_closure = f(D, A_seed; R_closure)) [BR2-01, BR2-02]

- One listed item (effect / closure structure) is now downstream of D + seed + rule instead of independently supplied.
  That is a bookkeeping compression **inside** A_res, not the elimination of a residual component.
- Time orientation is still not selected, and Lorentz structure is still not derived.
- Empirical status: **zero confirmed distinctive GRUT quantitative predictions** (B5 not reached).

## 8. Is B2 now sharpened enough to run?

**Yes.** B3 supplies the two firewalls B2 needs:
1. **The value / structure firewall (B3-6).** State-dependent MI, redundancy and distinguishability are **values on a
   fixed access structure**. They are never to be counted as access compression, so H_corr cannot be silently booked as A.
2. **The access variables are now explicit:** seed, partition, resolution, time. B3-9/10 shows record counts move with
   all four, so **B2 must hold A fixed** while varying the preparation.

**Proposed B2 question.** With A fixed, is GRUT's preparation-relative integrated arrow the same price as SCOUT's
H_corr|Σ (RENAMING)? The arrow in question is S6 = NET-ARROW-CONFIRMED (conservative parent + declared preparation).
Alternatively, is the price carried by the GRUT environment / noise layer (RELOCATION, cf. B1-4)? Does any GRUT structure
fix the preparation?

**B2 guardrails:**
- the D0 scope repairs (E-03 / E-04: continuous Rényi and continuous MI only);
- O-6 stays falsified, not reopened;
- no ħ / outcome / gravity.

## Status

- **BRIDGE REPAIR 01:** applied (BR1-01 … 03). **B3:** DONE.
- **TRUE COMPRESSION:** 0. **CONDITIONAL COMPRESSION:** 3 (+1 identity-grade). **RELOCATION:** 3.
- **CANONICAL UPDATE CANDIDATES:** B1-CUC-1, B3-CUC-1 (neither applied).
- **HARD STOP. Awaiting owner review before B2.**

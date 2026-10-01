# SCOUT_0 FRONTIER QUEUE — live routes (updated after each probe)

Not ranked by total score. Status / strongest objection / next decisive test per route.

| ID | Route | Status | InfoGain | DepReach | EmpReach | AssumpCost | ImportRisk | Strongest objection | Next decisive test |
|----|-------|--------|----------|----------|----------|------------|------------|---------------------|--------------------|
| P-06 | CM/non-CM boundary, Lorentzian-exponent family | **COMPLETE — success (a)**: S5-1 scaled kernel CM; terminal fully-CM member; p=2 proven counterexample; bracket [1.0,2.0] | high (serves S-7) | medium (generator layer only) | none | very low | low | bracket not yet a theorem (small-t only for 1<p<2) | P-06b proof of small-t mechanism; P-06c second deformation direction |
| P-08 | faithful-representation lift audit | QUEUED | high | high (lift layer) | none | very low | low | "full observable algebra" may be unstated for the classical substrate | read L0_LIFT_SELECTION_01.md survivor definitions; state faithfulness precisely; test |
| P-09 | locality-of-representation constraint | QUEUED (after P-08) | high | high | none | low (NEW ASSUMPTION declared) | low | locality of a lift may be ill-defined without a declared geometry on the quantum side | freeze definition; test survivors |
| P-17 | Caldeira–Leggett bath as C-B candidate | QUEUED | high | high (noise-origin seam) | none | medium | medium | may just reproduce standard QBM (= S-1 equivalence class, the designed null) | compute reduced mean response at O(t³); compare to C-B identity |
| P-18 | ergodic fast-variable null | QUEUED | medium-high | medium | none | low | low | averaging may kill noise entirely (null outcome is still informative) | theorem-grade statement on martingale survival |
| P-01 | floor-edge portability (Laplacian families) | QUEUED | medium | medium (substrate layer) | none | very low | low | floor edges may be definitional to C1-a and trivially port or trivially fail | spectral argument on symmetric vs asymmetric Laplacians |
| P-02 | boson-statistics sector formation | QUEUED | high | high (SFG-0 wall) | none | low-medium | low | free boson IR may be dominated by BEC-like zero mode, muddying classes | compute IR classes for two scaling sectors |
| P-23 | decoherence-rate scaling exponent | QUEUED (after P-06/P-17) | high | medium | **HIGH** | medium | medium | hidden supplied scales (geometry, localization, preparation) masquerading as "parameter-free" exponent | symbolic scaling with declared input inventory; DP/CSL comparison table |
| P-15 | nonlinear branching → Born weights | QUEUED | high | high (outcome frontier) | none | medium | medium | any mismatch with Born kills the model, not helps GRUT (by charter) | worked toy model with computable branch weights |

## Wave-2 spawn queue (from Wave-1 results so far)

- **P-06b:** prove the t^α small-t CM-killing mechanism for the full kernel → closes [1.0,2.0] as a theorem.
- **P-06c:** two-scale Lorentzian mixtures — is the Lorentzian exponent the *general* d_mono order parameter or family-specific?

# BRIDGE-1 PROVENANCE

| source | location | commit | role |
|---|---|---|---|
| **SCOUT reviewed** | `ryangrvr/TestingGRUT` branch `scout-2-reviewed` | `1f08ae156a346f6089ae3b3aab2be80b9ee4ef8f` | reviewed residual (`SCOUT_2_REVIEW_HANDOFF.md`, `review/REVIEW_ACCOUNTING_LEDGER.md`) |
| **Canonical GRUT** | `ryangrvr/GRUT-RAI` branch `master-w25bu9` | `b935099f61008acf111a762a6e346d947c2d0c38` | READ-ONLY comparison source |

## Canonical files consulted (read-only; fetched as a detached checkout with push disabled)

| file | what it pins |
|---|---|
| `uploads/GRUT_Working_Theory_2026-09-30.md` (primary) | supplied ledger, earned ledger, negatives, minimal architecture, claim fence |
| `EA0_ENDOGENOUS_ACCESS_FORMULATION_01.md` (with the `EA0_CORRECTIONS_01.md` precedence banner) | access candidates C-0 … C-9, L1 – L8 |
| `EA0_OWNER_RULING_02.md` | EA-0 = UNFORMULABLE at earned Level-0 scope. **§3: the local net is admissible substrate structure, existence / uniqueness not derived. §4: "H is k-local" = CRITERION-SMUGGLED** |
| `L0_1A_CHARTER_01.md` | the C1-a substrate: chain, N = 24, K = 0.3·I + L (unit springs), retained site 0 + bath |
| `L0_1B_CHARTER_01.md` | weighted networks w_ij = c_α\|i−j\|^{−α}, K = 0.3·I + L(w); the geometry predicate (resistance distance range Q ≥ N/2 and line ordering) |
| `L0_1C_CHARTER_01.md` | V = ½xᵀK_b x + βΣxᵢ⁴ (the on-site convex quartic) |
| `L0_ACCESS_BRIDGE_01.md`, `S2_NOISE_ORIGIN_HARD_DDET_01.md` | L0-1e noise Q = 2·diag(T_i) |
| `GS1_GEOMETRY_SELECTION_VERDICT_01.md` | geometry SELECTED-IN-CLASS from full site-resolved access; single-site access UNDERDETERMINED |
| **B3 additions:** `P5_ACCESS_VERDICT_01.md` | access = representational + closure-given-seed + access-split; the seed is declarative |
| `P6_SEED_SELECTION_VERDICT_01.md` | seed = NON-DERIVED INPUT; selected only up to block decomposition; L-A / L-B/E / L-C / L-D1 / L-D2 |
| `L0_ACCESS_BRIDGE_OWNER_RULING_02.md` | NONUNIQUE-LIFT; direct classical branch TRIVIAL/IDENTITY (end-site readout globally observable); scope §5 |
| `calc/p5_access.py`, `calc/p6_seed_selection.py` (**inspected, not run or imported**) | the closure *rules*: P-5 uses alg{I, H, B}; P-6 uses the unital product span of the seed, ad_H-closed, with `adh=False` for the L-D2 sufficiency test |
| `GRUT_WORKING_THEORY_SYNTHESIS_01.md` (Ledger A-13), `SYN0_OWNER_RULING_01.md` §2 | how access is booked: "access seed / declared readout / access sets" (A-13), layer 3 = "access / readout" |

## Precedence and use rules

- Later owner rulings govern over historical wording (e.g. the EA-0 corrections over the EA-0 formulation text).
- No canonical code is imported. The canonical models are **re-implemented** in `bridge/b1/` and `bridge/b3/` from the charter definitions.
- Claims are not copied between the two sources without a crosswalk classification.

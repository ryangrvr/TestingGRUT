# RESIDUAL SYNTHESIS STATUS

> **GRUT RESIDUAL SYNTHESIS 01 — FROZEN CONSOLIDATION OF EXISTING RESULTS**
>
> This is not a scientific campaign. It claims **no new physics result**.

## Pins

| item | pin |
|---|---|
| Branch base | `gravity-scout-1-frozen @ df25f29bf217931b366c0ef676cb38b7b34ad115` |
| Frozen inputs | `scout-2-reviewed @ 1f08ae1…`; `grut-bridge-1-frozen @ 753b90a…`; `qft-scout-1-frozen @ 5b07573…`; `gravity-scout-1-frozen @ df25f29…` |
| Canonical GRUT | `GRUT-RAI/master-w25bu9 @ b935099f61008acf111a762a6e346d947c2d0c38` (read-only, untouched) |
| **Freeze commit** | **the commit that introduces this file** on `grut-residual-synthesis-1` (parent `df25f29bf217931b366c0ef676cb38b7b34ad115`). A commit cannot contain its own hash, so the exact hash is pinned by **`grut-residual-synthesis-1-frozen`** |

## Queue — CLOSED

| item | file | status |
|---|---|---|
| master residual ledger | `MASTER_RESIDUAL_LEDGER.md` | DONE |
| constraint / coupling map | `CONSTRAINT_COUPLING_MAP.md` | DONE |
| supplied-input table | `SUPPLIED_INPUT_MASTER_TABLE.md` | DONE |
| information-flow ledger | `INFORMATION_FLOW_LEDGER.md` | DONE |
| selector challenge | `SELECTOR_CHALLENGE.md` | DONE (entry conditions only; **no selector campaign opened**) |
| master synthesis | `GRUT_RESIDUAL_SYNTHESIS_01.md` | DONE |
| integrity check | this file | DONE |

## Internal consistency checks performed

- Every final grade in the master ledger matches its source campaign's frozen / corrected wording. No grade is stronger
  than the frozen campaigns earned.
- **TRUE COMPRESSION = 0** across Bridge-1, QFT-SCOUT-1 and GRAVITY-SCOUT-1. SCOUT's scoped TRUE DERIVATIONS are kept in
  their own ledger and are not counted as residual compression.
- **Key statement:** supported only in its qualified form. Every uniqueness result is conditional on a supplied
  criterion / class / state / seed (`GRUT_RESIDUAL_SYNTHESIS_01.md` §12).
- **Toy results** (G2 ε ladder, charge-sharp products, G3 time-band conditioning) appear only as illustration-grade edges.
- **The heuristic d\*** appears only as a HEU edge, conditional on L_all.

## Integrity

- No frozen scientific branch changed: scout-0, scout-1, scout-2, scout-2-review, scout-2-reviewed, grut-bridge-1,
  grut-bridge-1-frozen, qft-scout-1, qft-scout-1-frozen, gravity-scout-1, gravity-scout-1-frozen.
- Canonical GRUT unchanged.
- **No PR. No merge. No selector branch opened.**

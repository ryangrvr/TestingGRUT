# SCOUT-2 INDEPENDENT REVIEW (post-freeze)

- **Branch:** `scout-2-review`, created from the frozen SCOUT-2 head `af0042f`.
- **The frozen `scout-2` ref is never modified.** Nothing under `results/`, `zoomouts/`, `ledgers/`, `probes/` or
  `SCOUT_2_HANDOFF.md` is edited on this branch.
- **New material lives only under `review/`.** Proposed wording repairs to frozen files are recorded as **errata**
  (`review/ERRATA_PROPOSED.md`). They are applied only if the owner unfreezes.

**Independence caveat.**
- The theorem review (IR-01 … IR-04) is the **owner's** independent hostile review.
- The numerical reproductions here were written by the same agent that ran SCOUT-2. They are **independent code
  paths**, not an independent reviewer:
  - fresh implementations;
  - different numerical routines;
  - fresh random frames.

  They can catch implementation errors and frame-sampling luck. They cannot catch shared conceptual errors in the
  diagnostics.

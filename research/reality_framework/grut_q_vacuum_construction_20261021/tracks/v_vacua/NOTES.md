# Track `v_vacua`: notes (index)

21-22 October 2026. This track is the independent verifier of `../vacua/`. The record is **`VERIFY.md`**:

| Section | Content |
|---|---|
| §0 | summary |
| §1 | verdict table, 29 rows |
| §2 | every correction with its corrected value |
| §3 | q-sector derivations |
| §4 | nucleation and SM numerics |
| §5 | logic of the P_q readings |
| §6 | credit |
| §7 | inputs and provenance |
| §8 | limits |
| §9 | reproduction |

`results/verdicts.json` holds the table in machine-readable form.

Labels are as in BRIEF.md: [D] [C] [Q] [E] [A] [R]. Every comparison with measured SM parameters is retrospective.
Downloads (SMDR v1.3, arXiv abstracts) are in the scratchpad (`.../qtheory/v_vacua_dl/`), and none is in the
repository. The track's PDG/TeX downloads in `.../qtheory/vacua_dl/` were read only.

**Key outcome.** No item is REFUTED. 20 are CONFIRMED and 9 CORRECTED.

The main correction:
- The stability/criticality line moves from **170.94 to 171.03 GeV** at (M_h, α_s) = (125.20, 0.1180), with SMDR
  matching and running. This agrees with Hiller et al. 2024 (171.04).
- The cause is the 2013 matching: y_t lacks the 4-loop QCD term (−0.00094), and g3(M_t) sits 0.0017 above the 4-loop
  QCD value.
- The retrospective z values become 3.15 (direct M_t), 2.68 (with the 0.3 GeV MC-to-pole error) and 1.73
  (cross-section mass), against the track's 3.34, 2.84 and 1.84.
- The held-out FCC-ee reach becomes 11.6σ (track 12.2σ).

The track's own pipeline reproduces to 2 MeV. Its β functions match SMDR's at 1-3 loops.

Reproduce: `cd code && sh v3_build_smdr.sh && python3 -B run_all.py` (VERIFY.md §9).

# Main-session checks of the report rows drawn from ledger entries s10–s12

20 October 2026. The full adversarial re-reads of s10–s12 did not complete in the first passes
(session limits); s10's later completed (its corrected entry is in `ledger/`). Before that, the
main session checked the specific claims that WORKING_THEORY.md's s10–s12 rows use, by grep
against the cited sources. This file records what was run and what it found, including what it
missed.

## What was run

```
# s10 (M1)
grep -n '0.987510\|0.935065\|0.998458\|0.999680\|778468\|778483' \
  research/reality_framework/grut_memory_candidate_20261019/{REPORT,OPEN_CLAIMS,WORKED_REGIME}.md
grep -n -i 'repetition code\|transverse-field Ising\|Peierls' \
  research/reality_framework/grut_memory_candidate_20261019_review/REVIEW.md
grep -n -i 'every sampled\|201 sampled\|all 201' research/reality_framework/grut_memory_candidate_20261019/{REPORT,OPEN_CLAIMS}.md
# s11 (toric, unmerged branch, read with git show)
git show origin/claude/toric4d-memory-vs-scrambling-20261019:research/reality_framework/toric4d_memory_vs_scrambling_20261019/REPORT.md \
  | grep -n '977\|1.64e+04\|7.98e+04\|0.77\|1.18\|post-hoc\|does not bear\|Independent review\|estimator'
# s12 (program ledgers)
grep -n -i 'FORBIDDEN\|+16\|empty' NO_GO_LEDGER.md POSTULATE_MAP.md
grep -n -i 'excluded by an input\|+16 net\|derived. tier' GRUT_FINAL_SYNTHESIS_01.md
grep -n -i 'nothing selects the access seed\|currently irreducible' GRUT_PROGRAM_STATE_SYNTHESIS_01.md
grep -n 'UNFORMULABLE' EA0_OWNER_RULING_02.md; grep -n 'NONUNIQUE-LIFT' L0_ACCESS_BRIDGE_01.md
python3 -c "import json; print(len(json.load(open('provenance/claims.json'))))"   # 74
```

## What it confirmed and corrected

- s10: the retention, readout and writer numbers match `REPORT.md:25,30,88-91` and
  `WORKED_REGIME.md:78,84`. Corrected: "retains better at every sampled time" to "no worse at
  every sampled time" (`REPORT.md:28`); the repetition-code reading attributed to the
  independent review (`REVIEW.md:36`).
- s11: the half-life numbers, the post-hoc face-spin correction, the estimator change, "does not
  bear on B", "Independent review outstanding" and the Metropolis/Davies-type bath match
  `REPORT.md:3,21,23,41-43,58-59,106,128` and `README.md:25,58`. Added: the fallback reading on
  O1 + O3 if the post-hoc correction is rejected (`REPORT.md:128`).
- s12: the 74-node register (53 + 21, `POSTULATE_MAP.md:6-7`), the absence of any banked
  FORBIDDEN no-go (`NO_GO_LEDGER.md:13,79`), "excluded by an input, never by a principle"
  (`GRUT_FINAL_SYNTHESIS_01.md:171`), the empty derived tier at +16 (`:236`), "currently
  irreducible" (`GRUT_PROGRAM_STATE_SYNTHESIS_01.md:73`), UNFORMULABLE and NONUNIQUE-LIFT. Corrected:
  "zero banked no-gos" to "no structural-impossibility (FORBIDDEN) no-go".

## What it missed (caught by the round-2 attack, now corrected in the packet)

- s10: the writer's energy flows from the register into the carrier, not the other way
  (`REPORT.md:116-118`); the original wording came from `REVIEW.md:42`.
- s12: the +16 count and the 53/21 split are the 2026-08-17 register sync stamps of the
  responsive-vacuum stance, not a current price of every model instance; "currently
  irreducible" predates the EA-0 ruling by a day; EA-0 is UNFORMULABLE because the earned
  substrate is classical, and the four-way candidate pattern is a finding over ten examined
  quantum candidates, not a no-go.
- s10: "Radiative stability is not established" was dropped at 4e9a161 and is restored in the
  packet's not-established list.

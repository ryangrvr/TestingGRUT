# Handoff-compliance audit of packet `grut_closed_whole_theory_20261019`

Lens: HANDOFF-COMPLIANCE AUDITOR. Read in full: `NEXT_CHAT_HANDOFF.md`, `START_NEXT_CHAT.txt`,
README, WORKING_THEORY, FORMAL_CORE, CALCULATIONS, OPEN_CLAIMS, SOURCES, MANIFEST, the
CAMPAIGN_LEDGER summary table and s01 entry, ledger s01/s10/s11 in full, targeted sections of
s03/s05/s09, `code/verify.py`, `code/finalize.py`, `results/m1_matched_control_rerun.txt`.
Executed: `python3 code/verify.py` and `python3 code/verify.py --rerun` (both exit 0, worst
relative deviation 0.0; log at `.../wf3/round1/handoff/verify_rerun.log`). Nothing under
`/home/user/TestingGRUT` was modified. Packet paths below are relative to
`/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019/`;
`handoff` = `research/reality_framework/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md`.

Repository state at audit time: branch `codex/grut-closure-fork-20261019`, HEAD `1489773`
("WIP checkpoint 7"), eight WIP commits after the handoff commit `96a7e7e`, all pushed
(`git fetch` then `git rev-list --left-right --count origin/...HEAD` = `0 0`).

## Findings (most severe first)

### F1 — BLOCKING — The posted report is a self-declared draft that forbids citation
- File / location: `WORKING_THEORY.md:1-2`:
  "DRAFT CHECKPOINT (19 Oct 2026, in progress). Sections marked `[[...]]` await the
  campaign-ledger workflow ... Do not cite this file until the banner is removed."
  Contradicted by `README.md:30-32` ("Status in one line: a working closed-system framework
  ... no additional law ... has survived").
- Problem: the handoff's deliverable is "the working theory", "not another handoff or a plan to
  start later" (`handoff:22-25`; `START_NEXT_CHAT.txt:3`). What is on the authorized branch at
  origin is a file that tells readers not to cite it. The banner is also stale: no `[[...]]`
  section remains inside WORKING_THEORY.md (grep), while the promised work it refers to is
  recorded elsewhere as pending (F2, F3).
- Deciding check: `grep -n '\[\[' WORKING_THEORY.md` -> only line 1;
  `git log -1 --format=%h origin/codex/grut-closure-fork-20261019` -> `1489773` (banner present).
- Severity: blocking (for the claim that the deliverable is finished and posted).
- Minimal fix: remove the banner only after F2–F4 are resolved; until then do not present the
  packet, in README or the coordination log, as the delivered working theory.

### F2 — BLOCKING — The ledger's own text contradicts the packet's claim that every sector was adversarially re-read
- File / location: `CAMPAIGN_LEDGER.md:3-7` (header: "Each entry below was written by one reader
  ... and then re-read adversarially by a second"); `WORKING_THEORY.md:188` ("Every entry was
  written ... and re-read adversarially by a second"); `WORKING_THEORY.md:345` (same, §8);
  `OPEN_CLAIMS.md:97` (same, verification table).
- Problem: `CAMPAIGN_LEDGER.md` prints, for eight of twelve sectors, "Verification: adversarial
  re-read PENDING (reader entry only)" — lines 780 (s05), 890 (s06), 1018 (s07), 1131 (s08),
  1246 (s09), 1369 (s10), 1480 (s11), 1609 (s12). Only s01–s04 carry correction counts
  (lines 176, 349, 516, 645). The packet therefore asserts a verification standard for the
  whole ledger that its own ledger file says was met for a third of it.
- Deciding checks:
  - `grep -n '^Verification:' CAMPAIGN_LEDGER.md` (12 lines; 8 say PENDING).
  - `scratchpad/wf1_combined.json`: `verified_ok = None` for s05–s12, `False` with
    `n_corrections` 9/12/8/11 for s01–s04.
  - `cmp ledger/ledger_s0N_*.md scratchpad/wf1/ledger_s0N_*.md` -> identical for all twelve
    (s05–s12 are the unrevised reader drafts).
  - Commit message `7d5c0bc`: "s01-s04 adversarially re-read; s05-s12 re-read pending".
- Severity: blocking (an exact-claim-boundary failure on the packet's own process claims,
  which the owner explicitly values; it also means the s05–s12 numbers in WORKING_THEORY §5
  rows and §6 rest on single-reader entries).
- Minimal fix: either execute the s05–s12 re-reads and rebuild, or change
  `CAMPAIGN_LEDGER.md:3-7`, `WORKING_THEORY.md:188,345` and `OPEN_CLAIMS.md:97` to
  "s01–s04 re-read adversarially; s05–s12 single-reader entries, not yet re-read".

### F3 — MAJOR — Four unfilled template placeholders are promises, not results
- File / location: `OPEN_CLAIMS.md:97` `[[LEDGERVERIFY]]`, `:98` `[[CRITIC]]`, `:99` `[[ATTACK]]`;
  `SOURCES.md:33` "Primary texts actually inspected for this packet (beyond the repository):
  `[[LITERATURE]]`".
- Problem: these rows sit in a table titled "Verification performed for this packet" and in the
  sources list, and are hashed into `MANIFEST.json` as final content. The workflow record shows
  the underlying work did not run: `scratchpad/wf1_result.json` has `literature: []`,
  `critic: None`, `formal: None`.
- Deciding check: `grep -rn '\[\[' *.md` -> the four placeholders (plus one legitimate `[[...]]`
  math bracket at `ledger/ledger_s07_wm_b1_b2.md:235`).
- Severity: major.
- Minimal fix: replace each with the actual record, or delete the row and state "not performed"
  in OPEN_CLAIMS §Verification (the handoff asks that verification be logged as actually
  performed, `handoff:41-42`).

### F4 — MAJOR — Primary-literature theorems are applied without the inspection the handoff required
- File / location: `FORMAL_CORE.md:216-221` (COND 5.3, Fewster–Verch natural-state exclusion
  with stated hypotheses "dynamical locality with non-trivial relative Cauchy evolution");
  `FORMAL_CORE.md:62-69` (DEF 1.5, Pollock et al. process tensor); `FORMAL_CORE.md:122-136`
  (THM 2.5, Feynman–Vernon); `FORMAL_CORE.md:250` (§7 uses COND 5.3 as the cost of any
  state-selection rule); `WORKING_THEORY.md:314-318` (same use).
- Problem: `handoff:197-202` says of exactly these three papers "Only search/abstract metadata
  was checked for this synthesis; inspect the relevant primary text before applying its
  theorems." The packet applies them and records no inspection (`SOURCES.md:33` is the
  `[[LITERATURE]]` placeholder; `wf1_result.json` literature list is empty).
- Deciding check: `SOURCES.md:33`; `handoff:197-202`; `python3 -I -c "import json;
  print(json.load(open('scratchpad/wf1_result.json'))['literature'])"` -> `[]`.
- Severity: major (uncertain whether the statements are wrong — COND 5.3 matches the
  Fewster–Verch abstract as I recall it — but the handoff instruction is unmet and the packet
  cannot show the hypotheses were checked).
- Minimal fix: inspect the three primary texts, record paper/section/theorem and hypotheses in
  SOURCES.md, and adjust COND 5.3's hypotheses if needed; or label COND 5.3, DEF 1.5 and
  THM 2.5 "cited from abstract metadata; hypotheses not checked against the primary text" and
  weaken §7's "must survive COND 5.3" accordingly.

### F5 — MAJOR — Independent-verification claims have no artifact in the posted record
- File / location: `OPEN_CLAIMS.md:92-93` ("Independent re-execution by a separate agent |
  5030 + 24 numeric leaves bit-exact"; "Independent re-implementation ... code unseen | 16
  quantities agree"); `CALCULATIONS.md:215-240` (same, plus two "adversarial code review"
  verdicts); `WORKING_THEORY.md:340-343`; `OPEN_CLAIMS.md:28-30` (Established item 6).
- Problem: nothing in `results/`, `code/` or the repository records these runs; a reader of the
  branch cannot check them. The claims are not fabricated — the artifacts exist outside the
  repository (`scratchpad/wf2/rerun/REPORT.md`, `wf2/reimpl/REPORT.md`,
  `wf2/review1/REVIEW.md`, `wf2/review2/REVIEW.md`) — but the handoff's standard is to "log
  verification actually performed" (`handoff:41-42`), and the posted packet offers only the
  assertion. The reproducibility claim itself is true: `code/verify.py --rerun` executed here
  passes with worst relative deviation 0.0 on both JSON files.
- Deciding check: `ls results/` -> three files, none a verification report;
  `git grep -l 'reimpl' HEAD -- research/reality_framework/grut_closed_whole_theory_20261019`
  -> none; scratch reports at the paths above.
- Severity: major.
- Minimal fix: add the re-execution and re-implementation reports (verdict sections plus the
  comparison JSON) under `results/verification/`, re-run `finalize.py`, and cite them; or
  reword to "re-executed in this session (verify.py --rerun)" and drop the separate-agent
  claims.

### F6 — MAJOR — Coordination-log completion entry not written; the draft contains more placeholders
- File / location: `AGENT_COORDINATION.md` (repo root). `handoff:38-42`: "Read
  AGENT_COORDINATION.md before editing, claim the intended area ... and log verification
  actually performed."
- Problem: `git diff 96a7e7e HEAD -- AGENT_COORDINATION.md` is exactly one added line, the
  CLAIMS row (line 29). No completion entry exists for the packet, although eight packet commits
  are pushed. The prepared entry (`scratchpad/log_entry_draft.md`) still carries
  `[[VERIFY_SUMMARY]]`, `[[ATTACK_SUMMARY]]`, `[[FINAL_COMMIT]]` and asserts "Each ledger entry
  was ... re-read adversarially by a second", which F2 shows is false as of HEAD.
- Deciding check: `git diff 96a7e7e HEAD --stat -- AGENT_COORDINATION.md` -> `1 insertion(+)`.
- Severity: major (an explicit handoff instruction; harmless to fix, but the draft must not be
  appended as written).
- Minimal fix: after F2–F5, append an entry stating exactly which re-reads, reruns and reviews
  ran, with their artifact paths, and the final commit hash; remove the claim row when done.

### F7 — MINOR — Count and wording errors in the "closure forces nothing local" claim boundary
- File / location: `OPEN_CLAIMS.md:15-17`: "Fourteen 'closure forces X' claims are each
  refuted by a legal closed model or a standard theorem (FORMAL_CORE §4); ten of them are
  exhibited with numbers in calculation 1 (S2–S7)". `FORMAL_CORE.md:180` titles §4 "(each with
  a countermodel)".
- Problem: the §4 table has 14 rows; nine, not ten, cite calc 1, and they span S1–S7 (row 1
  cites S1). Rows at `FORMAL_CORE.md:193-196` (Lyapunov exponent, compact topology, constants,
  cosmological energy conservation) give scope statements rather than countermodels; the
  Lyapunov row cites neither a model nor a theorem, so OPEN_CLAIMS' "each refuted by a legal
  closed model or a standard theorem" is not exact for it.
- Deciding check: `sed -n '/^## 4\./,/^## 5\./p' FORMAL_CORE.md | grep -c 'calc 1'` -> 9.
- Severity: minor.
- Minimal fix: "nine of them (S1–S7)"; retitle §4 "(each with a countermodel or a stated
  reason)"; for the Lyapunov row cite THM 2.9 (almost periodicity in finite dimension) or mark
  it scope-only.

### F8 — MINOR — WORKING_THEORY promotes the O-6 ruling to a demonstration the ledger says was never run
- File / location: `WORKING_THEORY.md:195` (s03 row): "the one closed-system item (O-6,
  conservative chain with internal bath) shows closure plus partition plus product start does
  not force a strict arrow (D-4 recurrence obstruction)".
- Problem: `OPEN_CLAIMS.md:66,78` and `ledger/ledger_s03_level0_floor.md:101-105,168,191`
  record that the O-6 chartered run on the declared chain was never executed ("Nothing was
  computed on the declared pin-0.3 chain"; "not frozen, not run, and banks nothing"); the
  FALSIFIED ruling rests on identities, the D-4 theorem and a non-member chain (pin 0.45,
  N = 800), and `ledger_s03:214` lists promoting those non-member numbers as declared-chain
  results under "must not be promoted". "Shows" in the row does exactly that compression.
- Deciding check: `ledger/ledger_s03_level0_floor.md:191` and `:214`; `OPEN_CLAIMS.md:78`.
- Severity: minor.
- Minimal fix: "the O-6 ruling (identities, the D-4 recurrence theorem and a non-member chain;
  the chartered declared-chain run was never executed) ...".

### F9 — MINOR — Two sentences state V7's content as fact rather than as the handoff's account
- File / location: `WORKING_THEORY.md:104-108` ("it is where GRUT V7's parent object actually
  lives"); `FORMAL_CORE.md:134-135` ("This places the program's V7 'genotype' object:
  downstream, conditional, and standard").
- Problem: `ledger/ledger_s01_v7_ctp.md:41` (assumption 6) and `:125-126` (discrepancy 2) say
  the identification of V7's S_CTP with the record's S_IF "is the handoff's ... not the
  record's; it cannot be verified here" because `GRUT_V7.md` is at no path in any of 992
  commits. Elsewhere the packet is careful (`SOURCES.md:30`, `OPEN_CLAIMS.md:70`,
  `WORKING_THEORY.md:182`).
- Deciding check: `git log --all --name-only | grep -i v7` -> only
  `GRUT_III/calc/v7_access_geometry_memory.py`.
- Severity: minor.
- Minimal fix: "as the handoff describes V7's parent object (the file is not in the
  repository)".

### F10 — MINOR — The toric-bridge refutation is listed among withdrawn additions without its required qualifiers
- File / location: `WORKING_THEORY.md:47-49` (§0: "the self-correction/no-scrambling bridge in
  the classical 4D toric code" in the list of additions that "did not survive");
  `WORKING_THEORY.md:296-298` (§6 bullet).
- Problem: the result comes from an unmerged branch whose code was not run here
  (`ledger/ledger_s11_toric_bridge.md:9`), whose report says "Independent review outstanding"
  (`tb:REPORT.md:3`, cited at `ledger_s11:112`), and whose O2 verdict depends on a labelled
  post-hoc observable correction (`ledger_s11:138-140`). §7 and §8 carry "not merged" and "no
  toric code run", but §0 and §6 present it on the same footing as results executed in the
  repository. The handoff: "read before importing ... Preserve that distinction"
  (`handoff:124,194-195`).
- Deciding check: `ledger/ledger_s11_toric_bridge.md:9,112,138-140,173-183`.
- Severity: minor (the ledger entry itself is exemplary; only the summary lines lack the label).
- Minimal fix: append "(as reported on the unmerged branch; not re-run here; independent review
  outstanding)" at both places.

## Confirmed (handoff items met, with the check)

- All handoff campaign rows covered. The handoff table has ten rows (`handoff:113-124`:
  V7/CTP, influence hierarchy, GRUT I–II, GRUT III, CP-12, WM/B1/B2, Candidate 1,
  closure/no-boundary, M1, toric), not nine; they map to s01, s02, s04–s11. s03 (Level-0) and
  s12 (program no-go/postulates) are additional.
- Every `ledger/ledger_s*.md` has question, assumptions/model class, executed evidence with
  line citations, answer, closure verdict, classification, remaining freedom, failures
  retained, must-not-promote, handoff-row discrepancies, unread/missing, sources (heading grep
  over all twelve).
- Corrected M1 commit used: `WORKING_THEORY.md:6-7`, `SOURCES.md:9`, `ledger_s10` (pinned to
  4e9a161/96a7e7e/da8f925, byte-identical). All handoff numbers (`handoff:138-162`) match
  `results/m1_matched_control_rerun.txt` and `ledger_s10 §3b-3f`; "matter readout 0.778468"
  in the s10 row is M1's value (`grut_memory_candidate_20261019/REPORT.md:89`), the control's
  is 0.778483. Facilitation benefit stated as withdrawn (§0, §6, §7, OPEN_CLAIMS).
- Candidate 1 kept conditional (`WORKING_THEORY.md:200` s08 row; §7 "for its supplied
  couplings"; OPEN_CLAIMS conditional list); tree and loop orders kept separate; no
  determinant recomputed (§8).
- C2: "adds no selection condition beyond restriction, trace matching, positivity,
  normalization and modular KMS" (§6); B three-prescription control "registered, unexecuted"
  (§6, §7, OPEN_CLAIMS, `ledger_s09 §3h`).
- No cross-import: `WORKING_THEORY.md:302-309` states Candidate 1, M1 and GRUT III live in
  different state spaces with no common action/state space/constraints/calibration.
- No new handoff file; no new branch; the eight commits stage only packet files plus the one
  coordination claim row (`git show --stat` of each); branch pushed and in sync with origin.
- Toric report treated as read-not-run: `ledger_s11:9` ("No branch code was run"; two result
  JSONs parsed with stdin-only `python3 -I`), `SOURCES.md:28`, `WORKING_THEORY.md:347`.
- V7 absent from the repository and declared so: `SOURCES.md:30`, `OPEN_CLAIMS.md:70`,
  `WORKING_THEORY.md:182`, `ledger_s01 §1`. No private file IDs anywhere in the packet (grep).
- Not called validated universal physics; gravity/"all reality" not announced derived:
  `README.md:30-32`, `WORKING_THEORY.md:60-63`, `OPEN_CLAIMS.md:34-47`.
- Completion standard: §0 answers "what GRUT now says / what calculates / what closure
  constrains / which additions did not survive / what remains supplied" in five labelled
  paragraphs.
- Core deliverables present: explicit theory datum (`FORMAL_CORE.md` DEF 1.7; CALCULATIONS
  calc-1 table); CTP downstream (THM 2.5); preparation-to-measurement with energy/correlation
  accounting (S1); non-monotone relaxation (S2, S5 exact return at T = 4π); non-unique local
  law (S5); initial-correlation example with the inhomogeneity ¼Z₀ in closed form (S4);
  nested/overlap compatibility with the Bell-pair chain (C1, C2) and quantitative edges
  (C3, C4); predictions conditional on inputs (§7); the extra condition a distinct prediction
  would need (§7; FORMAL_CORE §7 table with costs and the known-physics baseline).
- Reproducibility: `code/verify.py` (hashes + 32 checks) and `--rerun` both pass here; worst
  relative deviation 0.0 on 5030 + 24 numeric leaves.
- Spot-checked program numbers quoted in WORKING_THEORY: net +16 (`POSTULATE_MAP.md:6-7`);
  ISW 32σ -> 2.0σ (`GRUT_II_What_Survived.md:75`); Γ_T 62.7 orders (`ledger_s04:85`, verifier
  run); 0.338 coherence split (`ledger_s02:85`, P3_NC_LIFT_VERDICT); Ω_K = −1/(e^{2N}−1)
  arithmetic (−4.54·10⁻⁵, −2.06·10⁻⁹); DESI DR2 3.1σ is a retained record
  (`SIGNATURE_AUDIT.md:56` via `ledger_s04:146`), not a resumed cosmology dispute.
- No new postulate proposed; FORMAL_CORE §7 gives the cost/baseline table the handoff requires
  for any that might be.

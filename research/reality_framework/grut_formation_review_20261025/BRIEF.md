# Brief: independent review of Codex's formation, receiver and conversion packets (issue #13)

25 October 2026. Claude's independent review was requested in Codex's CLAUDE_HANDOFF.md files of the four packets
below. Read the packets only through git (`git show <sha>:<path>`, or `git archive <sha> <dir> | tar -x -C
<scratch dir>`). Never edit or check out Codex's branches, and never write inside Codex's packet directories.

| Commit | Branch | Packet | What it claims |
|---|---|---|---|
| bb211d1 | origin/codex/grut-core-conversion-20261023 | research/reality_framework/grut_core_conversion_20261023/ | Hidden conversion channel X̄ + core(Q) → X + core(Q−2); radial convergence (35.15% retained) |
| 59a1a71 | origin/codex/grut-conversion-streaming-20261023 | research/reality_framework/grut_conversion_streaming_20261023/ | First-conversion histories; ω/m_X = 0.861195; horizons 195 pc / 2.37 Mpc / 8.38 Mpc; mass fractions 3.175%, 3.359%, 6.718%, 5.469% |
| 9636950 | origin/codex/grut-soft-receiver-20261023 | research/reality_framework/grut_soft_receiver_20261023/ | The old boundary reflects about 32% at p/m_X = 0.02; the antiparticle estimate is absorber-sensitive (29%) |
| ff49e5d | origin/codex/grut-causal-receiver-20261024 | research/reality_framework/grut_causal_receiver_20261024/ | Causal CFS receiver: reflection < 2e-8 (charged), 6.6e-4 (neutral); retention 35.36%; extraction closes to < 1e-6; nine checks, 48 arrays |

Read each packet's CLAUDE_HANDOFF.md (ff49e5d and 59a1a71 list the exact items to challenge), REPORT.md,
MATHEMATICS.md and code/.

**Context.** This is Codex's empirical branch of GRUT: a Higgs-portal-extended working model whose charged scalar
objects (Q-ball-like cores) form from a supplied condensate. Candidate 1 is unchanged and paused. The vacuum
density is underived. Everything here is one model with supplied fields; nothing is evidence that the fields exist.

## Rules for every agent

- **Where to write.** Write only inside `research/reality_framework/grut_formation_review_20261025/tracks/<key>/`:
  NOTES.md (method, numbers, verdicts, limits, reproduction command), code/ (your own code), results/*.json. Large
  scratch output goes in /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/<key>/.
- **Git.** Read-only only (show, log, archive). No commits.
- **Independence.** Independent means your own code and method. Running Codex's code to reproduce its outputs is a
  separate, labelled check (reproduction). Say which one each number comes from.
- **Verdict per checked claim.** One of: CONFIRMED / CORRECTED (give the corrected value) / REFUTED / NOT CHECKED. Give
  your number next to Codex's.
- **Claim labels.** Derived [D], computed [C], quoted [Q], estimate [E], argued [A], recalled [R].
- **Scope.** Never claim validation, exclusion or a new theory. Comparisons with existing data are retrospective.
- **Credit.** Credit published methods: CFS/PML absorbers (Berenger; Kuzuoglu–Mittra; Roden–Gedney), Q-ball
  formation and fragmentation (Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki), and mixed or decaying dark
  matter structure constraints.
- **Return value.** A digest of at most about 700 words: verdict table, key numbers, corrections, and what remains
  uncertified.

# Brief: what does closure force that standard physics does not? (parallel pass, 27 October 2026)

Owner request (issue #13 context): "what does the closure principle force that standard physics doesn't?" This pass
runs in parallel with Codex's dark-matter work and must not touch Codex's packets.

## Where the program stands (read these first; do not repeat them)

- `research/reality_framework/grut_closed_whole_theory_20261019/FORMAL_CORE.md` and `OPEN_CLAIMS.md`.
  - Closure (DEF 1.1: a closed whole with no exterior; every subsystem's bath is the rest) forces nothing beyond
    standard quantum theory in the finite class.
  - Of 14 "closure forces X" claims, 6 are refuted by countermodels, 4 are true by definition and 3 are out of scope.
- `research/reality_framework/grut_vacuum_principle_20261021/SYNTHESIS.md` (the vacuum pass).
  - The obstruction: physics is invariant under H → H + cI on a fixed background.
  - Every closure-motivated route to the vacuum energy is a relabelled published proposal: q-theory, sequestering,
    unimodular, everpresent Λ, Banks–Bousso and others.
  - Kaloper–Padilla collapse was judged not testable. DESI DR2 tests were already done there.
- `research/reality_framework/grut_q_vacuum_construction_20261021/HANDOFF.md` (draft, **unreviewed**).
  - Klinkhamer–Volovik q-theory with the final-state equilibrium condition P_q.
  - Under the natural readings, P_q requires the electroweak vacuum to be **absolutely stable**, a Standard Model
    relation on m_t(m_h, α_s). Measured masses lie on the wrong side at about 3.1σ (retrospective).
  - The relation agent never ran, so the prediction was never frozen.
- Codex's foundation gate (git show 2908e91:research/reality_framework/grut_concept_checkpoint_20261025/FOUNDATION_GATE.md).
  - The null result is retained.
  - Compactness, finite state capacity, an observer dictionary and vacuum equilibration are **extra assumptions**,
    not consequences of the word "closed".
  - Any proposal must state:
    1. the exact additional postulate;
    2. compatibility;
    3. a baseline counterexample it forbids (a strict restriction S_C ⊊ S_base, with a derivation);
    4. one fixed discriminator and comparator, chosen before tuning;
    5. the supplied inputs and a failure condition.
- Codex's exact architecture (git show 04c53f5:research/reality_framework/grut_exact_architecture_20261026/THEORY.md and
  REPORT.md).
  - SM + metric + φ, X, Y, Σ + hidden U(1) (charges 1, 3, 6). A compact T³ branch is a stated choice.
  - Includes a closed clock-and-record model and the scalar operator basis.
- The GRUT dark-matter branch:
  - research/reality_framework/grut_formation_review_20261025/REVIEW_*.md;
  - the heavy-mediator spot point: H_i = 1.3335×10¹¹ GeV, M_Y = 30 H_i, h = 40 H_i (Σ vev scale 5.33×10¹² GeV),
    κ = 3.4587×10⁻³, e = 10⁻⁷;
  - Candidate 1 low-energy target: m_X = g v, g = 0.1, v = 10 m_X, λ_X = 10⁻³, λ_φ = 0.03, m_X ≈ 92 MeV;
  - Higgs-portal obligations in REVIEW_REHEATING.md;
  - our one-loop running code: research/reality_framework/grut_formation_review_20261025/tracks/hw_gates/code/gate4_running.py.

## The task

Take candidate **strengthenings** of closure: extra postulates that closure motivates. Put each through the
foundation gate, and quantify it against current data and against the GRUT branches. Find out whether any of them
gives GRUT distinctive, testable content. If none does, say so.

For each postulate, report:
- the exact postulate;
- what it forces, with a derivation;
- the baseline it forbids;
- the discriminator and comparator, with numbers and sources (label retrospective comparisons);
- literature credit, and whether it is GRUT-native or a relabelling;
- its consequences for the GRUT branches: Candidate 1, the heavy window, the dark-matter chain;
- a verdict: DISTINCTIVE-TESTABLE / STANDARD / ALREADY EXCLUDED / UNTESTABLE / INCONSISTENT WITH GRUT BRANCH;
- the failure condition.

## Rules

- Write only in `research/reality_framework/grut_closure_forcing_20261027/tracks/<key>/`: NOTES.md, code/,
  results/*.json, nothing over 2 MB. Put scratch in
  /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/closure/<key>/.
- Git is read-only. No commits.
- Use claim labels [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled. Cite primary
  sources with arXiv numbers. For current data use WebSearch or WebFetch (load them with ToolSearch "select:WebSearch,WebFetch"),
  and quote values with the source.
- Never claim a validation, an exclusion of the real world or a new theory. Mark comparisons with existing data as
  retrospective. A postulate adopted by choice is not a consequence of "closed".
- Return a digest of at most about 800 words.

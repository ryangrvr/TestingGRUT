# Brief: review of Codex b73de64 (conditional growth versus source clocks), issue #13

27 October 2026. Codex's packet `research/reality_framework/grut_source_timing_20261026/` at commit b73de64 (branch
origin/codex/..., read through git only) extends 42 scalar/longitudinal phase-tangent histories to τ = H_i t = 800 on
the two quiet reheating backgrounds (m_I, Γ_I)/H_i = (3, 0.5) and (5, 1). It claims:

1. A supplied 0.03 rad preparation-phase seed reaches unit fractional |X| modulation by sampled τ = 134 (a = 23.053)
   and 147.5 (a = 20.878). A 10⁻⁶ rad seed does so by τ = 772 and 782. A 10⁻⁹ seed does not before τ = 800. The local
   charge modulation reaches 1 earlier (τ 120, 134.5 for 0.03 rad).
2. Fixed late widths give decay-number medians τ = 6.89×10⁶ (Y_L), 4.36×10¹⁷ (vector, e = 10⁻⁶) and 4.36×10²⁰
   (e = 10⁻⁷). The late radiation constant C = a²H/H_i is 1.986 and 1.476. At the 0.03 rad unit samples, the Y_L
   rest-clock decayed fraction is 1.35–1.48×10⁻⁵ and its injected conformal-energy fraction is 3.72–4.31×10⁻⁸.
3. Therefore, for the supplied seeds, linear control is lost long before the prescribed decay heating. The next
   calculation must carry an early nonlinear release into the heating era.
4. The auxiliary quartic orbit has μ/k ≈ 0.440 as k → 0, a band edge at 3.135–3.140 and a second band at 5.57–5.98.
   The wavelength-cube charge at k = 2.13 is 1.34–2.09×10⁵ (geometric only).

Codex's request (CLAUDE_HANDOFF.md):
- start with MATHEMATICS.md and verify_response.py;
- check the direct-vector adaptation against the full transfer, the phase-current sign convention and the matched
  radiation clocks;
- check whether the conditional time ordering is stated correctly, without reading the linear extension as a
  nonlinear history.

The rules of [BRIEF.md](BRIEF.md) apply (write only in `tracks/<key>/`; git read-only; no commits; label independent
code versus reproduction; a verdict per claim; claim labels [D] [C] [Q] [E] [A] [R]; never claim validation,
exclusion or new theory; return a digest of at most about 700 words).

## Prepared environment (already done; do not redo)

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
# Codex tree (git archive of b73de64, the whole research/reality_framework directory):
#   $SP/b73/research/reality_framework/
# Pinned kinetic inputs that Codex's ks_common bootstrap needs (from e71730c):
#   $SP/pinned/research/reality_framework/grut_formation_review_20261025
export GRUT_KINETIC_INPUT_ROOT=$SP/pinned/research/reality_framework/grut_formation_review_20261025
export GRUT_KINETIC_WORKDIR=$SP/kwork GRUT_ONSET_WORKDIR=$SP/owork
export GRUT_TIMING_WORKDIR=$SP/review/<key>/twork GRUT_MODE_WORKDIR=$SP/review/<key>/mwork   # your own
cd $SP/b73/research/reality_framework/grut_source_timing_20261026/code && python3 -B -c "import timing_common"
```

Only the track `st_rerun` may run Codex scripts that write into `$SP/b73` (they save JSON into the packet directory).
It must first copy the published JSON files it will overwrite to `$SP/review/st_rerun/published/`. Every other track
treats `$SP/b73` as read-only: import from it or copy files out of it.

Useful own material (Claude, already committed):
- `tracks/rp_reheating/code/rp.py`: our own reheating engine, which reproduced Codex's e160356 runs. Its `extra`
  group has yr_b3_g0.5 and yr_b5_g1 in `tracks/rp_reheating/results/extra.json`, with yields 1.0454206 and 1.0474062.
- `tracks/rm_review/code/quartic_floquet.py`: our own quartic-orbit Floquet code (orbit rebuilt from E and L).
- `tracks/kin_source/results/source_interface.json`, `tracks/hw_gates/results/gate1_vector.json`: our source widths and
  release X̄ bracket (0.13–2.6 X̄ per net charge).
- Reviews: REVIEW_REHEATING.md, REVIEW_REHEATING_MODES.md, REVIEW_COREWAVE.md.

Machine: 4 cores, shared by up to four tracks. Use at most 2 threads or processes each, and keep runs under about
60 minutes.

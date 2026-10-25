# Track `critic`: adversarial consolidation

Claude review, 25 October 2026. The deliverable is [CONSOLIDATED.md](CONSOLIDATED.md), which holds the verdict
table, the cross-track disagreements, the overreach list, the corrections for Codex, what is certified and what is
open, and the next calculation.

## Method

- **Inputs.** I read every `tracks/*/NOTES.md`. Five are consolidated: reproduce, independent_solver, spectrum,
  operators and conversion. xcharge covers packet 8986d34, which is outside this brief. affleck_dine has no notes.
- **Codex packets.** Read through `git archive` into
  `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/critic/`:
  - CLAUDE_HANDOFF.md of ff49e5d, 59a1a71, bb211d1 and 9636950;
  - the ff49e5d REPORT and run JSONs;
  - bb211d1's RATE_RESULTS.json.
- **Own analysis.** `code/critic_checks.py`. I ran no Codex code and no new PDE evolutions. The inputs are:
  - the raw histories of the independent_solver track (its own spectral solver);
  - the spectrum and conversion result files;
  - Codex's JSON.
- **Where the numbers come from.** Every recomputed number is [C] from that script. The [Q] inputs are named in the
  script.

## Recomputed numbers [C]

| Check | Result |
|---|---|
| K1: old damper minus box, Q60 (pp of Q0) | +0.037 at τ = 2000, −0.588 at 2500, −0.309 at 4000, −0.280 at 9600, −0.059 at 12000; mean over 9600–12000 −0.188 |
| K1: a 1e-6 amplitude perturbation | window-mean shifts ≤ 0.024 pp; max \|ΔQ60\| ≤ 0.027 pp |
| K2: escape through r = 240 by τ = 12000 under that perturbation | ΔQ = +13.6, ΔE = +11.7 |
| K3: five converged runs from three codes (Codex P0125/P00625, spectrum F2dt/F3h, independent_solver h = 0.5, dt = 1/64) | retention 35.3565–35.3631%; Q_esc 23515.2–23526.9 (mean 23519.2); E_esc 25478.8 and 25486.1 |
| K3: Codex I025 against K3 | Q +0.69%, E +0.54% |
| K4: late-window offset correction [E] | 34.574 + 0.161 = 34.735%; independent_solver 34.740%; spectrum 34.741% |
| D: N_X̄ over five converged values | 255.0–256.1, mean 255.5 (0.656% of Q0) |
| K5: rate table, PCHIP (log) against the dense scan | ×1.431 over p = 0.45–0.56 (×1.643 with linear PCHIP); ×1.097 over the full band, which equals the optical-depth ratio 54.61/49.74 = 1.098 |
| K6: free_soft packet (\|u\|² ∝ exp(−(r−300)²/100²), k0 = 0.02) | energy fraction 0.300 below k = 0.0163 and 0.079 below k = 0.01; the return-coverage threshold estimate is k ≈ 0.016 [E] |

**Return times [A].** I take a wave emitted at r = 60 and reflected near r ≈ 300 in the old 270–360 layer. It comes
back after about 680 time units at p = 1, 1070 at p = 0.5, 4800 at p = 0.1 and 9600 at p = 0.05.
- So the deficit that turns on at τ ≈ 1800–2500 must come from fast waves.
- The recovery inside the window fits p ≈ 0.05.

## Limits

- No new formation runs were made. The mechanism analysis rests on one pair of runs, independent_solver's box run
  and its old-damper run.
- The late-window and τ = 48000 corrections assume constant offsets [E].
- The structure-screen critique is methodological. I did not rerun CLASS.

## Reproduce

    python3 -I research/reality_framework/grut_formation_review_20261025/tracks/critic/code/critic_checks.py \
        research/reality_framework/grut_formation_review_20261025/tracks/critic/results/critic_checks.json

The paths to the scratchpad runs and the Codex extractions are set at the top of the script.

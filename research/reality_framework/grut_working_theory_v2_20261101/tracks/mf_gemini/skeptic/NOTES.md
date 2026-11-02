# mf_gemini skeptic: working notes

The write-up is `../SKEPTIC.md`. This folder holds the skeptic's own code (`code/`) and outputs (`results/`). None of it
imports from `../code`. Same-session AI work; every shape, coupling, criterion and convention is a model choice; no
validation, exclusion or GRUT prediction is claimed; data comparisons are retrospective.

| script | what it recomputes | output |
|---|---|---|
| `code/sk_g1g2.py` | G1/G2 by numerical solve; loop caps from renormalisation conditions; Higgs routes (both signs, Landau pole, Higgs loop on χ); large-N_f and the species bound; the −v domain light-Y mass; K_crit with the Y deposits rescaled | `results/sk_g1g2.json` |
| `code/sk_background.py` | LSODA background; M1r regression K = 1.985805; slices (S1, ε_H = 1, v-crossing, rest point); K for sudden and constant-Γ decay in three normalisations | `results/sk_background.json` |
| `code/sk_waterfall.py` | independent linear modes (e-fold time, adiabatic start, v = a^{3/2}u) + nonlinear δN; σ dropped at gχ = 100H; band-exit detection; 12 cases in quad, lin, hill and log shapes | `results/sk_waterfall.json` |
| `code/sk_ns.py` | direct N_*; hilltop window; consistent finite-σ_c log family with exact Mukhanov–Sasaki n_s and α_s | `results/sk_ns.json` |
| `code/sk_pbh.py` | f_PBH by entropy conservation; slow-reheating formation; no-waterfall patches | `results/sk_pbh.json` |
| `code/sk_preheat.py` | KLS crossing map with a consistent φ decay into radiation | `results/sk_preheat.json` |

Run with `python3 -B` from `code/`. Per-case waterfall files go to scratch
(`$SP/wt2/mf_gemini_skeptic/wf/`) before `combine`.

Self-found bugs, both fixed before any reported number:
1. Mode starts placed deep inside the horizon (fixed: start at max(k/aH = 50, background start), heavy adiabatic
   vacuum).
2. u-variable underflow during long heavy pre-histories, which made explicit RK45 unstable (fixed: integrate
   v = a^{3/2}u; all cases re-run).

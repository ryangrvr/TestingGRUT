# Review handoff: lattice QED in GRUT detector histories

From: Claude (Perplexity Computer), branch `claude/grut-qed-detector-histories-20261107`.
Base: `origin/codex/grut-photon-bridge-review-20261107` at `e88f3ba`.
Owner-authorized: work through the repository by branch and PR, and hand off to ChatGPT/Codex.

## What this packet adds

The packet extends the photon bridge (`PHOTON_MATCHING.md`, `PHOTON_BRIDGE_REVIEW.md`) from free photons to interacting charged matter. The model is the lattice Schwinger model (1+1D QED) with Kogut–Susskind fermions and the Gauss law solved exactly. It tests the integrated paper's matching condition (14) outside free fields. Everything here is written independently of the GRUT codebase and imports nothing from it.

| File | Content |
|---|---|
| `ed/PREREGISTRATION.json` (+ `.sha256`) | Original predictions P-A through P-E, written before any computation |
| `ed/schwinger.py`, `ed/spectrum*.py`, `ed/extrapolate.py` | Exact diagonalization, N ≤ 24 |
| `ed/histories.py` → `ed/histories_results.json` | GRUT two-time Weyl-instrument histories, compared with the two-point (quasi-free) prediction |
| `mps/PREREGISTRATION_AMENDMENT_01.json` (+ `.sha256`) | Sealed in `dd1e107` before any production DMRG. Targets and windows are unchanged. |
| `mps/schwinger_mps.py` | Exact Gauss-law MPO, bond dimension 6 |
| `mps/validation_vs_ed.json` | Gate: matches ED to rel ≤ 2e-13 on 4 cases |
| `mps/analyse.py` | Recipe from the amendment, committed in `c4e1ddf` before results were read |
| `mps/production_m*.json`, `mps/DMRG_RESULT.json` | DMRG grid and verdicts |
| `lightcone.py` → `lightcone_results.json` | Wilson-fermion vs plaquette-photon \(O(k^2a^2)\) velocities; LHAASO Crab electron ceiling |
| `particles.py`, `gravity.py` (+ results) | Earlier probe: Nielsen–Ninomiya doubling, Wilson doublers, lattice Newton law, spin-2 light bending |
| `GRUT_Gravity_Particles_Probe.md`, `GRUT_Play_It_Out.md` | Plain-language reports |

## Results

- **P-C, P-D, P-E pass (ED).** GRUT histories of the smeared electric field approach the two-point prediction at m=0 (deviation 0.0076 → 0.0041 as the lattice is refined). They stay away from it at m/g=0.5 (about 0.009, not shrinking). The first instrument measurably disturbs the second.
- **P-A and P-B miss under ED** (+3.3% and +9.9%). These remain on record and are not overwritten.
- **P-A and P-B pass under the sealed DMRG rerun.** P-A gives 0.5641906 against 0.5641896 exact. P-B (holdout) gives 0.519232 against 0.51918. All 12 fit variants lie within 0.25% of their targets, and the chi 160/120 change is below 5e-8. See `RESULTS.md`.
- **Light cone.** Wilson electrons are faster than plaquette photons in every direction sampled. Crab electrons at ≥ 2.3 PeV then give \(a \lesssim 3\times10^{-32}\) m.

## Suggested bounded review (for Codex)

1. Rerun `mps/validate_vs_ed.py`, then at least one production point. Report what you reran, not what you read.
2. Check the MPO channel algebra in `schwinger_mps.py` against `ed/schwinger.py` by hand. The convention is sz = 2Sz, q_n = Sz_n + (−1)^n/2, and the hopping sign is +x.
3. Check that `analyse.py` implements the amendment's recipe exactly. Any deviation voids the verdict and must be logged as such.
4. Check the P-C/P-D caveat. Var(cF) differs between masses (0.24 vs 0.12). A matched-variance rerun is the clean follow-up.
5. Do not widen windows or change the recipe and call it the same test. A new recipe needs a new sealed amendment.

## Open, not attempted here

- **Chiral obligation.** Which escape from Nielsen–Ninomiya GRUT adopts (overlap/GW, domain wall, mirror decoupling or symmetric mass generation), and whether access or closure selects it.
- **Universal massless spin-2 coupling** as a derivation target. GR1 already marks it as an import.
- **A shared emergent light cone** across photon, fermion and any graviton, without hand tuning.

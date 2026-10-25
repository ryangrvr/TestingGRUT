# Track sp_critic: notes (consolidation of sp_point, sp_yledger, sp_vector, sp_phaseprep for Codex 373d37e)

25 October 2026 · Claude (issue #13). The verdicts, ledger and recommendations are in [CONSOLIDATED.md](CONSOLIDATED.md);
this file records method, numbers, limits and reproduction. The Codex packet was read only through
`git archive 373d37e research/reality_framework/grut_inflation_scan_20261025 | tar -x` into the scratchpad. No Codex
code was run and nothing was committed. It is one supplied tree-level model: not evidence, validation or exclusion.

## Method

1. I read the four tracks' NOTES.md and results JSON, Codex's CLAUDE_HANDOFF, REPORT, MATHEMATICS, PROTOCOL,
   CORRECTIONS, FINDINGS and Y_SOURCE_RESULTS, and our prior REVIEW_GAUGEPOINT and gp_critic CONSOLIDATED.
2. I cross-checked every shared number between tracks and against the packet.
3. I recomputed the disputed or load-bearing numbers with my own code in `code/`:
   - `critic_checks.py` [D, C]:
     - A: Y-decay energy under the three conventions;
     - B: the combined ledger, assembled from the track JSONs with the pair-to-f conversion re-derived
       (allowance = 0.01 B/(2 m_X Y_Q) = 0.0049296, f = 2.0285 × pairs, cold remainder 1.0143 (1/ε − 1));
     - C: an independent longitudinal-Proca mode ODE testing the frozen superhorizon gauge angle;
     - D: the decay-averaged Y energy and Y pairs for all 2,203 published Codex grid releases (Codex columns [Q]),
       compared with the benchmark thresholds and with Codex's own 1e-9 regularity test;
     - E: the ensemble-averaged Y energy of Codex's 31 published robustness ensembles, compared with the case-L
       threshold.
   - `rmax_relic.py` [C]: our sp_yledger/gp_leftovers maximal-equilibrium benchmark, scanned in injected debris energy
     R = E_inj/E_X, giving the f = 0.01 thresholds for case H (homogeneous) and case L (lumps plus conversion).
   - `b1_ppi_slope.py` [C inputs, E conversion]:
     - the noise-averaged response d⟨θ_Σ(∞)⟩/dψ0 under phase noise of 0.026 rad (the residual spread of
       sp_phaseprep's c_S = 1 preparation), with 5 offsets × 256 common draws, computed with sp_vector's vrel.c;
     - its conversion to B1 vectors through sp_vector's P0 transfer, with the preparation's blue phase spectrum.

## Key numbers

All of these are in `results/*.json`.

- **Y energy at decay.** 274.85 (H = Γ), 388.7 (t = 1/Γ) and **348.99** (decay-averaged, N_Y,eig, m_L) times E_X.
- **Frozen gauge angle.** The ratio of P1 to P0 is 2.356e-6, against (m_vac/m_inf)² = 2.3607e-6, for k = 3e-6 to
  3e-5. This confirms sp_vector piece A.
- **Ledger.**
  - Debris: 0.131–0.140 pairs per net charge, f = 0.27–0.28, which is 26.5–28.3 times the allowance.
  - Frozen total: f = 3.47–3.49.
  - Benchmark equilibrium: f = 15.2–15.9 (sp_yledger).
- **Case H thresholds.** f = 0.01 at R = 0.014, 0.69, 2.06 and 7.3 E_X for ε = 0.19, 0.2402, 0.3267 and 0.598.
- **Case L lump thresholds.** f = 0.01 is crossed between R = 10–30, 30–50, 75–100, 100–150 and 150–200 E_X for
  ε = 0.0677, 0.10, 0.15, 0.19 and 0.2402. At ε = 0.3267 it is not crossed up to R = 200.
- **Scan.**
  - 161 of the 2,203 releases have \|ε\| ≥ 0.19.
  - Their Y-only R_Y is 16–2506 E_X (median 454), and their Y pairs are 0.0041–2.4.
  - Only 3 are regular (R_Y = 47, 161, 344), and only 1 is robust under σ = 0.03.
  - All 161 fail case H; the smallest ratio of R_Y to the threshold is 7.8.
  - 16 are nominally below the case-L threshold, but only one of these is regular (8.03e12 GeV).
- **Codex's 31 robustness ensembles.** None is below the case-L threshold. The closest are candidates 10 and 17 (c_H = 10,
  4.35e11 GeV), at R_ens = 115 and 124 against 80–109 and 72–97, counting Y alone.
- **B1 under the c_S = 1 preparation.**
  - Slope: d⟨θ_Σ(∞)⟩/dψ0 = 3624 ± 310 rad/rad.
  - Size at k_*: δψ(k_*) = 1.16e-3 rad per ln k, so δθ_G(k_*) = 0.70 against the criterion 0.226.
  - Result: 0.062–0.117 pairs per net charge, which is 13–24 times the allowance (8.7–32 across ±2 SE).

## Limits

- **Homogeneous k = 0 releases only.** No lattice was run, and fragmentation is expected near a ≈ 600.
- **The relic is a benchmark.** It assumes full thermalisation, 4→2 freeze-out, and case L with the bb211d1 conversion.
  It is not a kinetic calculation.
- **The thresholds assume the yield equals Y_Q at the light-scale onset (ymult = 1).** They come from a 6–9 point grid
  in R, so they are brackets.
- **Scan-wide R_Y omits the Σ and vector debris.** It is therefore a lower bound. It uses Codex's single-mass J_Y at
  a = 600–1000.
- **B1 uses phase-only noise.** The preparation's radial modes (1.40 H and 1.45 H) are left out. The P0 conversion
  assumes frozen superhorizon angles. The slope fit treats the offset points as independent, and a paired estimate is
  also reported.
- **The Proca check is a single-Higgs reduction.** It is exact for superhorizon modes.

## Reproduction

The commands are listed in CONSOLIDATED.md §9.

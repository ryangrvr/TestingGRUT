# Track `reproduce`: reproduction audit of ff49e5d (causal receiver), 9636950 (soft receiver) and the retention numbers of bb211d1

25 October 2026 · Claude Code review track · read-only on Codex's branches; no commits.

Two kinds of numbers appear below, always labelled:

- **REPRODUCTION** means Codex's own code, unmodified, run from a `git archive` extraction of
  ff49e5d (which contains the unchanged 9636950, bb211d1 and c4b99f3 packet trees). The extraction is in
  `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/reproduce/tree`.
  The pristine published copies are next to it in `ff49e5d/`, `9636950/` and `bb211d1/`.
- **INDEPENDENT** means my own code in `code/` and a different method from Codex's.
- **NEW CONTROL** means Codex's old executable driven by my own driver at parameters Codex did not run.

## 1. Environment, compiler and deviations from REPRODUCE.md

| Item | Codex (ENVIRONMENT.json) | This reproduction |
|---|---|---|
| Compiler | g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0 | g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| Flags | -O3 -std=c++17 (hard-coded in run.py, run_interface.py, refinement.py, operator_controls.py) | see below |
| Python / NumPy / SciPy | 3.12.14 / 2.3.5 / 1.17.0 | 3.13.16 / 2.5.3 / 1.18.1 (Matplotlib absent, so report.py and plot.py were not run) |
| CPU | not recorded | 4-core Intel Xeon 2.1 GHz, x86-64 baseline (SSE2, no FMA), shared with other review tracks |

**Flags [C].** The brief asks for -O2. `code/gxx_O2_wrapper.sh` sits first on PATH, rewrites Codex's `-O3` to
`-O2` and logs every compile line. -O2 builds were used for all of these:
- radial_nodal.cpp and radial4.cpp (bb211d1);
- operator_controls.cpp;
- the ff49e5d smoke run;
- the opening segments of P00625 (to τ=2810), P0125 (to τ=8040), P025_remote (to τ=2210) and P025_remote_time (to τ=1070).

For this code the -O2 build is about 3.4 times slower than -O3 (19.4 s against 5.65 s on the same 150-unit case).
The machine was also shared with another review track, so at 02:31 the wrapper was switched to keep Codex's
`-O3` for radial_cfs.cpp and radial_flux.cpp only. The cut-short -O2 runs were rerun from scratch at -O3.
The switch rests on byte-identical outputs from the two builds in three checks:
- smoke: all four raw-array hashes equal the published hashes;
- a 150-unit nonlinear run at h=0.25: final, history and memory CSVs are byte-identical;
- the -O2 opening segments are bitwise equal to the -O3 reruns: P00625 (282 history rows, 12035 trace rows), P0125
  (805 rows), P025_remote (222 rows) and P025_remote_time (108 rows).

This agrees with expectation: x86-64 has no FMA by default, ISO C++17 mode turns off contraction, and GCC never
reassociates floating-point sums without -ffast-math. The full compile log is in the scratch file `gxx.log`.

**No run length was reduced.** Every reproduced run used Codex's parameters and tmax.

**Reproduction gaps found [C].** These are not physics errors.
1. bb211d1's `refinement.py` never creates `data/`, and REPRODUCE.md does not say to create it. Both C++ programs
   write through unchecked `std::ofstream`. When the directory is missing, a run computes for its full length,
   writes nothing, and fails only when Python tries to read the files afterwards. I lost the first N00625/H003125
   attempt this way.
2. The nine old-spatial-method H runs of bb211d1 were **not rerun** (CPU budget). They are not inputs to any
   claim audited here. bb211d1's `convergence.py` needs all 18 runs, so its spatial and control statistics were
   recomputed from the nine N runs by `code/compare_reproduction.py` instead.

## 2. Reproduction result in one line

Every regenerated output equals Codex's published output **exactly**.

| Packet | What matched |
|---|---|
| ff49e5d | Every numeric leaf of all 18 regenerated JSON files equals the published value; only timing fields differ. The JSON files are the 12 run summaries, DISCRETE_LAYER, COLD_GRID, OPERATOR_CONTROLS, SOURCE_AUDIT, BOUND_MODE_DIAGNOSTICS and VERIFICATION. The last five plus DISCRETE_LAYER and COLD_GRID are byte-identical files. All 48 raw float64 array hashes match. |
| 9636950 | BOUNDARY_RESULTS, TRACE_AUDIT_RESULTS and VERIFICATION_RESULTS are byte-identical. |
| bb211d1 | The nine nodal summaries (N05 … N0125_far) match in every numeric field, and all nine trace SHA-256 values match. |

Details are in `results/reproduction_comparison.json`.

## 3. Verdict per claim

Labels follow the brief: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

| # | Claim (Codex value) | Reproduction [C] | Independent / new control [C] | Verdict |
|---|---|---|---|---|
| 1 | ff49e5d: discrete CFS reflection, charged band, below 2e-8 (max 1.82e-8 over 40 charged samples) | 1.8221975e-8, bitwise | Own operator with a source-driven lead fit: 1.8222e-8. All 16 rows with R>1e-14 agree to <1e-6 relative. The own operator matches the executable's exported C++ RHS to 2e-12. | CONFIRMED |
| 2 | Neutral (m=√6) maximum 6.6e-4 | 6.516076e-4, bitwise | 6.516076e-4 | CONFIRMED |
| 3 | Four free packets, error energy / initial energy: 1.1e-12, 4.32e-15, 3.78e-13, 6.88e-6; boundary-isolated (same-RK4) ≤ 5.5e-14 | identical | Snapshot coverage (own code): the reflected return can re-enter the measured region by tmax for only 80%, 80%, 93% and **70%** of packet energy. free_soft tests k ≳ 0.016 only. | CONFIRMED as numbers; scope CORRECTED (see §4.1) |
| 4 | 9636950: old damper reflects about 32% at p/m_X=0.02 (31.73%) | BOUNDARY_RESULTS byte-identical | Own transfer-matrix code: 31.735%. The other entries match too: 57.07, 16.71, 3.415, 0.1258, 4.921 and 14.83%. | CONFIRMED |
| 5 | Six formation runs: finite-time retention 35.36%; two finest grids within 0.019% | 35.3631102% / 35.3565086%; 0.01867%; observed order 4.68 | NEW CONTROL, old executable with its damper moved 270→720: 35.158%, **35.375%**, **35.383%** at h=.25, .125, .0625. At h=.125 with the damper at 1440: 35.335%. | CONFIRMED (35.36%); grid agreement understates the perturbation sensitivity (§4.2) |
| 6 | Time-step control (P025 vs P025_time) | +0.0305% retention; energy error 1.53e-5 vs 4.89e-7 (coarse step fails 1e-5, as disclosed) | — | CONFIRMED |
| 7 | Layer 360→720: retention shift 4.144e-7 relative; raw L2 inside 300: 1.25e-7 | 4.14407e-7; 1.24988e-7 | — | CONFIRMED |
| 8 | I025 has trace, final and memory bitwise identical to P025_time | identical | Physical history columns 0–7 are bitwise identical too (not claimed by Codex) | CONFIRMED |
| 9 | Extraction-cut closure below 1e-6 (energy 3.30e-7, charge 6.28e-8) | 3.2991896e-7 / 6.2813863e-8 | Analytic re-derivation of J_E and J_Q from the 7-point hopping matrix agrees with instrument.py [D] | CONFIRMED |
| 10 | 9636950: antiparticle estimate is absorber-sensitive, 28.95% | 28.95%, bitwise | Own estimators: 29% with Tukey windows, **+127% with a Hann window**, **×8–10 for a τ≥6000 window**. A net-flux estimator gives about 0% with Tukey windows and +32% with Hann. | CONFIRMED as computed; "29%" is not a stable measure of the sensitivity (§4.3) |
| 11 | bb211d1: radial convergence, 35.15% retained (N00625 35.15321%; spread 0.0693%; abs +0.231%; box +0.090%) | all bitwise | Moving the old damper to 720 raises the same executable's value to 35.383% (+0.654%). The old box control (R=360) moved it by only 0.09%. | CORRECTED: about 35.38% (old method, remote damper), which supports ff49e5d's 35.36% |
| 12 | Nine checks pass; 48 arrays verified | 9/9; 48/48 | These are author self-consistency checks, not independent tests | CONFIRMED (as reproduction) |

## 4. Corrections, nuances and undisclosed items

### 4.1 What the time-domain packets actually test [C]

Each packet error is a single snapshot at tmax. A Fourier component whose return from the wall or layer has not
re-entered r < R−4h by tmax contributes nothing to the error.

`code/packet_coverage.py` uses the lattice group velocity and the real CFS coordinate stretch inside the layer.
It finds that the round trip is complete by tmax only above these momenta:
- free_X: k ≥ 0.085;
- free_s: k ≥ 0.165;
- free_soft: k ≥ 0.0163. About 30% of free_soft's energy, including the 7.9% at k<0.01, cannot yet have reflected
  back.

So the "60000-unit low-momentum case" does **not** test p/m_X = 0.01–0.016 in the time domain. That end of the
band rests on the frequency-domain lead calculation alone, which I confirmed independently (1.82e-8 at k=0.01).
Codex correctly calls the packet numbers "not a pure reflection coefficient"; it does not state this coverage limit.

### 4.2 The finite-time retention is a perturbation-sensitive observable [C]

All runs of both discretizations agree to about 1e-5 in Q60(τ) until τ≈1400–1600. After that, differences grow
quickly:
- old damper at 270 vs 720: 2e-6 at τ=1400, 5e-4 at τ=2000, 7.6e-3 at τ=3000 (absolute charge fraction);
- grid differences grow in the same way.

The late-window mean reacts much less, but by more than 0.019%. Old-executable values at h=0.125 for different
damper radii R:

| R | 270 | 360 | 720 | 1440 |
|---|---:|---:|---:|---:|
| Retention | 35.145% | 35.177% | 35.375% | 35.335% |

These bracket the CFS value of 35.357%. The CFS layer itself is insensitive: moving it shifts retention by 4e-7.
The causal receiver therefore removes the boundary perturbation that the old method could not control.

Using the old-method results with a remote damper as an independent cross-family check, ff49e5d's 35.36% is
confirmed to within about 0.1%, not 0.02%:
- the two finest old/remote values (35.375, 35.383) lie 0.05–0.06% above the CFS pair (35.357, 35.363);
- the R=1440 value lies 0.06% below.

Fixed order-4 extrapolants are 35.3836% (old, remote damper) and 35.3636% (CFS). Neither is certified.

**Mechanism attribution [A].** Codex attributes the old 0.6% bias to "a boundary that later failed its soft-wave
test". Travel times do not support the soft-wave part. A p=0.05 wave (v≈0.05) needs about 12000 time units to go
from r=60 to the R=270 layer and back. The strongly reflected p ≤ 0.05 waves (3.4–57%) therefore cannot have
returned to the core within the run window. The divergence onset at τ≈1400–1600 instead matches returns of fast
waves, v ≳ 0.4 (p ≳ 0.4), where the old layer reflects a few percent (4.9% at p=0.5 and 14.8% at p=1). The bias is real and
boundary-caused; its source is mainly the old layer's high-momentum reflection.

### 4.3 The 29% sensitivity depends on the window and the estimator [C]

In 9636950, the conditional all-band antiparticle fraction for N0125 changes as follows:
- Tukey windows (α=.025–.1): 1.15% → 1.48% when the absorber is halved, i.e. +29%.
- Hann window: 0.448% → 1.016%, i.e. +127%.
- τ ≥ 6000 window: 0.038% → 0.404% (×10.6). With the same late window, the grid sequence N025/N0125/N00625 gives
  0.078/0.038/0.032%.

So the 1.1–1.5% fraction comes mainly from early-record content that the Tukey window keeps at full weight. A
net-flux estimator (no in/out split, own code) is absorber-insensitive under Tukey windows. Codex's conclusion,
that the number cannot supply an anti-X population, is reinforced. The specific "29%" should not be quoted as the
size of the instability.

### 4.4 Code-audit items not disclosed in the reports [C]

1. **Unchecked output streams.** Neither radial_cfs.cpp nor radial_nodal.cpp checks its `std::ofstream`.
   bb211d1's refinement.py does not create `data/`, so a fresh checkout loses a full run silently (§1).
2. **Hard-coded physics and observables in C++.** The coupling 0.1|X|² (λ term .025|X|⁴), the neutral quartic 3
   (m_s=√6), the portal S²|X|², the core radius 60 and the runaway threshold Fmax>100 are not parameters. All are
   consistent with the stated action (dV/dS, dV/dX and the energy density checked [D]). There is no silent
   clipping: the runaway check throws, and `clamp` only shapes σ(r).
3. **Late window.** The late window is τ>0.8·tmax, i.e. 9610–12000 (the text says 9600–12000). For P025_remote it
   is 19210–24000.
4. **Energy-ledger normalization.** Ledger errors are relative to the total E0≈38951, which is almost all rest
   mass (E0−mQ0 = −17.9). The absolute ledger error is 0.019. Compared with the escaping excess energy
   E_esc−mQ_esc ≈ 1938 at the I025 cut, that is 1e-5, still small. The relative errors are not errors in the
   radiation energy.
5. **The receiver acts on the neutral field too**, with the full quartic force left unstretched in the layer.
   The mass-shift monitor covers only the charged mass² (S²−1+0.1|X|²).
6. **verify.py thresholds are looser than the report's claims.** The cut closure check is <1e-5, the packet
   check <1e-3 and the grid check <1%. The claimed numbers do pass the tighter statements.
7. **The -O3 flags are disclosed** in ENVIRONMENT.json. -O2 gives byte-identical outputs on this platform.

## 5. What remains uncertified (unchanged by this audit)

- A positive asymptotic particle spectrum: the signed cut currents and the receiver work are not particle counts.
- A bound/coherent subtraction. The instantaneous bound-mode charge is resolution-dependent: 11183 at h=.25 and
  about 13320 at h=.125 and .0625.
- A nonlinear matching theorem for the layer, and reflection below k=0.01. The h=2 cold-grid value is already
  1.35e-4 at k=0.005.
- Asymptotic retention: Q60 is still decaying, 34.57% at the τ=19200–24000 window.
- 3D or cosmological formation, general hot-core leakage and the vacuum density.

None of this is evidence that the supplied fields exist.

## 6. Reproduce this audit

Scratch directory: `S=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/reproduce`.

```bash
# extraction (read-only git)
git archive ff49e5d research/reality_framework/grut_{causal_receiver_20261024,soft_receiver_20261023,core_conversion_20261023,daughter_core_20261023} | tar -x -C $S/tree
mkdir -p $S/tree/research/reality_framework/grut_core_conversion_20261023/data   # missing from REPRODUCE.md
# REPRODUCTION: Codex scripts with the -O2 wrapper on PATH (see code/run_queue.sh for the full job list)
code/run_queue.sh $S/tree $S
(cd $S/tree && python3 research/reality_framework/grut_soft_receiver_20261023/code/boundary.py)
(cd $S/tree && python3 research/reality_framework/grut_causal_receiver_20261024/code/{source_audit,bound_modes}.py)
(cd $S/tree && python3 research/reality_framework/grut_causal_receiver_20261024/code/verify.py --with-raw)
(cd $S/tree && python3 research/reality_framework/grut_soft_receiver_20261023/code/trace_audit.py && python3 research/reality_framework/grut_soft_receiver_20261023/code/verify.py --with-traces)
python3 -I code/compare_reproduction.py $S results/reproduction_comparison.json
# INDEPENDENT
python3 -I code/indep_old_damper.py results/indep_old_damper.json
python3 -I code/indep_cfs_reflection.py results/indep_cfs_reflection.json $S/build/operator_controls
python3 -I code/packet_coverage.py results/indep_packet_coverage.json
python3 -I code/indep_antiparticle.py $S/tree/research/reality_framework/grut_core_conversion_20261023/data results/indep_antiparticle.json
# NEW CONTROL (old executable, own driver)
python3 -I code/old_far_absorber.py $S/build/radial_nodal $S/far OLD0125_R720 0.125 0.03125 810 720   # likewise h=.25,.0625; R=1440 with L=1530
```

Results: `results/*.json`; compile lines: `results/compile_log.txt`; job log: `results/job_status_log.txt`.
Credits: CFS/PML absorbers (Berenger; Kuzuoglu–Mittra; Roden–Gedney); Q-ball formation (Kusenko–Shaposhnikov;
Enqvist–McDonald; Kasuya–Kawasaki).

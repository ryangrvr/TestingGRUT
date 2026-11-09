# B2-01r1 detector reconstruction: result

Completed 8 November 2026 America/Chicago (9 November UTC). **Geometry gate passed; overall massless matching remains unresolved because the lattice-to-matched gate failed.** No corrected massive point was evaluated.

The quadratic reconstruction was derived from cell-average and boundary constraints without fitting probabilities. It was frozen before production at commit d010ad185ac2fc99c42de5229b1a7653670905c0, following a documented import-only startup correction. The numerical production took 111.42 seconds, within the 2400-second cap. All four lattice sizes, the N14 refinement, seven timing protocols, mode doubling, baseline replay and post-run audit completed.

## What works and what remains unresolved

At N20, the maximum geometry TV over the original timings fell from 0.0010107531 (B1) to 0.00000653145 (B2), about 155 times smaller. Including the three new timing protocols, the maximum B2 geometry TV is 0.00000918199, below the 0.00025 gate. The new timings were specified before any B2 production history was inspected.

The largest lattice-to-matched TV is 0.00087507757 at [0,2], above the frozen 0.000740 allowance by 0.00013507757 (18.25%). The confirmation set also misses that component allowance: 0.00078021570 at [0,1.7]. Thirteen of fourteen gates pass, but every gate was required. Unused geometry allowance is not reassigned after seeing the result.

The direct target maximum is 0.00087689042 and the maximum additive total, including the 0.00001 numerical reserve, is 0.00088738427. Both are below 0.001. These finite observations are useful progress, but they do not change the failed frozen verdict or establish continuum convergence.

Removing a component error need not improve a previously cancelled total. On the original timings the maximum direct target TV increased from B1's 0.00042023649 to B2's 0.00087689042. At the first event it increased from 0.00000801062 to 0.00025925708. The new rule improves the independently measured detector geometry; it does not improve every lattice-to-target number.

## N20 probabilities summarized by total variation

| Times | B1 geometry | B2 geometry | B2 lattice-to-matched | B2 direct target | Additive + reserve |
|---|---:|---:|---:|---:|---:|
| 0 | 0.000241768331 | 2.35714663e-07 | 0.000259021362 | 0.000259257077 | 0.000269257077 |
| 0,1 | 0.000605546526 | 4.32751738e-06 | 0.00031234388 | 0.000308759634 | 0.000326671398 |
| 0,2 | 0.000540604164 | 2.30670309e-06 | 0.000875077565 | 0.000876890424 | 0.000887384268 |
| 0,1,2 | 0.0010107531 | 6.53145462e-06 | 0.000847788611 | 0.00084336986 | 0.000864320066 |
| 0,0.6 | 0.000519212434 | 3.74012278e-06 | 0.000259021362 | 0.000259257077 | 0.000272761485 |
| 0,1.7 | 0.000586457557 | 1.22533714e-06 | 0.000780215699 | 0.000779226076 | 0.000791441036 |
| 0,0.6,1.7 | 0.00102201406 | 9.18199151e-06 | 0.000690685137 | 0.000685677489 | 0.000709867129 |

## Finite-spacing trend and the next discriminating test

| N | Maximum B2 geometry TV | Maximum lattice-to-matched TV | Matched TV / a^2 |
|---|---:|---:|---:|
| 10 | 0.00014033374 | 0.00319485526 | 0.0154630995 |
| 14 | 4.65830386e-05 | 0.00172319432 | 0.0155087489 |
| 18 | 1.24199516e-05 | 0.00106856714 | 0.0154301095 |
| 20 | 9.18199151e-06 | 0.000875077565 | 0.0154363682 |

The ratio is nearly constant (0.01543-0.01551), consistent with an O(a^2) leading discrepancy on this finite grid. This is a post-run diagnostic, not a fitted correction or a proved asymptotic rate. The unresolved contribution can include Hamiltonian discretization, boundary matching and the assumed cell-average operator identification; this experiment does not separate them.

A focused next study would freeze the same detector and extend to N22, retaining these timings and thresholds. Simple a^2 scaling predicts approximately 0.00072951 at N22, only slightly below 0.000740; that is a forecast, not a passed result. A new protocol should budget the larger state space and an independent refinement, and predeclare what happens if the prediction fails. Do not change the weights to chase the observed error. N>20 was outside B2's resource card and was not run. A derivation or independent test of the lattice-to-cell-average dictionary remains necessary even if a finer point passes.

## Numerical and integrity checks

* Reproduction of constants, linear and quadratic profiles: maximum error 3.33e-16. Symmetry, weight sum and sine-coefficient quadrature checks pass.
* Synthetic oscillator ordered-Gaussian versus truncated Hilbert-space TV: 2.25e-16. Independently assembled N6 Hamiltonian and dense-versus-sparse full histories agree to 2.22e-16 and 1.49e-16, respectively.
* Fresh B1 replay: maximum TV 4.40e-15. N14 refinement changes B2 histories by at most 2.55e-14.
* Maximum eigensolver residual: 3.20e-11; probability normalization error: 2.63e-13; free-evolution norm error: 2.79e-13.
* Maximum 8192-to-16384 mode step: 1.72e-10; conservative analytic Gaussian history-tail bound evaluated in floating point: 5.60e-7.
* N18-to-N20 maximum TV: 0.00019492641, below 0.00025. Maximum two-event disturbance discrepancy: 0.00043565849, below 0.001.
* Post-run audit recomputed the assessment identically, checked every lattice/reference distribution and all stored ground-state residuals, and verified source hashes. Deliberately corrupted probabilities, a forged verdict, a changed source hash, and a simulated bypass of the massive prerequisite were rejected. The instrument precheck also rejects incorrect normalization.

The dense small-system checks are independently formulated; the full production and arithmetic audit reuse inherited code. Eigensolver residuals, refinement and a numerical reserve do not constitute interval-certified error bounds. These results concern only the supplied massless free-vacuum reference, detector, boundary convention and finite timing list. They neither validate all of GRUT nor establish an interacting-continuum or Standard Model result.

## Archive and reproduction

The source freeze, derivation, protocol, raw probabilities, decisions, audit script and checkpoint manifest are committed on the separate codex/grut-detector-reconstruction-20261108 branch. No merge, release or Zenodo deposit was made. Prior P1/E1/E2/M0/M1 verdicts are unchanged.

The delivered ZIP additionally contains all five binary ground-state checkpoints, original input/reference files needed by the audit, and SHA-256 hashes of its payload. Binary checkpoints are preserved in the ZIP rather than uploaded through the text connector. GitHub's CHECKPOINTS.json records their hashes and sizes. The source hashes in the archived run describe the exact Windows checkout bytes; the ZIP preserves those bytes, including line endings.

For a fresh production reproduction, check out the freeze commit in Git, install requirements.txt in an isolated Python environment, and run run.py with a new output directory outside the checkout. The runner requires a clean checkout and refuses existing output directories. For an offline audit of the delivered ZIP, install those requirements and run:

    python audit_b2.py --root source --output raw

The audit writes audit.json exclusively. Since the packet already contains the recorded audit, use a copy of raw with that one file removed when rerunning it; keep the supplied original archive intact. The machine-readable assessment.json is the authoritative frozen decision. A successful execution is compatible with a failed scientific gate.

# B2-01: frozen detector reconstruction study

Base: TestingGRUT commit 7d57fc7c8a36c8f349165afc946d9323f174e630 (B1-01, draft PR #14). This separate branch does not change P1/E1/E2/M0/M1 or reopen their certificates. The authoritative machine-readable parameters and thresholds are in protocol_b2_01.json.

## Before production

The exploratory B1-only diagnosis at N=20 compares four continuum functionals: true f; cell-averaged interior with true edge f; true interior with boundary-value edge replacement; both B1 approximations. Report both orders of the signed difference decomposition, with no percentage attribution. Only the four previously used timings are inspected in this stage. Raw output is exploratory/B1_geometry.json.

The single candidate is the quadratic moment reconstruction in DERIVATION.md. It has no fit parameters. Before freezing, verify degree-0/1/2 reproduction, partition, reflection, sine coefficients against independent quadrature, an oscillator reference, and an independently assembled N=6 Hamiltonian/history calculation at synthetic times 0,0.23,0.91. A wrong detector normalization must be rejected. These checks do not inspect the B2 production histories.

Commit this source, protocol, derivation, exploratory data and prechecks before running production. Run only from a clean checkout of that commit. Preserve the commit, source hashes, dependency versions, exact probabilities and ground-state checkpoints. One candidate and one production attempt are allowed; retain exceptions and resource-cap failures. A code correction requires a separately documented protocol revision, never silent overwriting.

## Registered runs and decisions

Use N=10,14,18,20 and an N=14 refinement. Keep L=5, sigma=0.8, detector coupling=2, phase=0 and m_R=0. Reproduce B1 using the same freshly computed states. Evaluate old times [0], [0,1], [0,2], [0,1,2] and new confirmation times [0,0.6], [0,1.7], [0,0.6,1.7]. The new timings are fixed before any B2 production probabilities are examined.

Use 8192 continuum modes and double to 16384. Refine N=14 eigensolver tolerance from 1e-11 to 1e-13 and split every evolution interval in two. Maximum per point: 600 seconds; total execution: 2400 seconds; one BLAS thread. No N>20 or parent-basis extension.

All gates must pass: prechecks; lattice and reference normalization/positivity; Gauss identity; eigensolver residual <=1e-9; B1 replay <=1e-9; refinement <=1e-8; reference mode doubling <=1e-7; analytic history tails <=1e-6. At N=20 every registered history must have geometry TV <=0.00025, lattice-to-matched TV <=0.00074, direct target TV <=0.001, and additive TV plus 0.00001 numerical reserve <=0.001. N18-to-N20 TV must be <=0.00025. Two-event disturbance discrepancies must be <=0.001.

Only if every massless gate passes, apply the unchanged detector at N=20, m_R=0.5, and compare old timings with existing P1 E14/q128, retaining every overflow outcome. Otherwise do not evaluate a corrected massive point. No interacting-continuum certificate is available from this transfer even if its finite comparison is small.

## Reproduction

Python 3.12, numpy==2.3.5, scipy==1.17.0. From a clean committed checkout, with outputs outside that checkout:

    python research/reality_framework/grut_detector_reconstruction_20261108/run.py --output /new/external/output/directory

The runner writes assessment.json, massive_status.json and execution.json, plus per-point raw JSON/NPZ/log files. A completed experiment may legitimately fail its scientific gates. A successful process exit does not mean the theory or continuum matching is established.

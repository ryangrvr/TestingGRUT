# DMRG rerun of P-A and P-B: both pass

The run follows `mps/PREREGISTRATION_AMENDMENT_01.json` (sha256 `0db11998…e18e`), committed in `dd1e107` before production DMRG. The analysis script was committed in `c4e1ddf` before the reported inspection of those results. Targets and ±3% windows are the same as in the original `ed/PREREGISTRATION.json`. This is a registered rerun of previously examined benchmarks, not a new untouched holdout. The original ED registration and results entered the repository together; its stated earlier local chronology is not independently established by that commit.

| Prediction | Target | DMRG primary | Error | Variant spread (12 fits) | Verdict |
|---|---:|---:|---:|---|---|
| P-A, vector mass at m=0 | 0.5641896 (exact, \(1/\sqrt\pi\)) | 0.5641906 | +1.7e-6 | 0.564179–0.564588 (target inside) | PASS |
| P-B rerun, binding energy at m/g=0.25 | 0.51918 (Byrnes et al. 2002) | 0.5192320 | +1.0e-4 | 0.519207–0.520345 (target 0.005% below lower edge) | PASS |

## How to read it

- **Fit-choice sensitivity.** The twelve-fit spread is 0.07% at m=0 and 0.2% at m/g=0.25. It is not a complete uncertainty bound. The agreement of the P-A primary to 2e-6 should not be quoted as the method's precision. Every tested variant lies within 0.25% of its target; untested finite-size or continuum systematics are not thereby bounded.
- **Bulk values at each lattice spacing** (primary finite-size fit, \(\omega-2m\)):

  | | x=9 | x=16 | x=36 | x=64 |
  |---|---:|---:|---:|---:|
  | m=0 | 0.6373 | 0.6192 | 0.6010 | 0.5919 |
  | m/g=0.25 | 0.5934 | 0.5754 | 0.5570 | 0.5477 |

  Both rows converge smoothly in \(1/\sqrt{x}\).
- **Why ED missed.** The N ≤ 24 lattices had physical lengths of only 4–14/g. The DMRG lattices are 12–24/g, up to N=192. The ED misses (+3.3%, +9.9%) stay on record as a failed first test. This rerun is a second, separately sealed test, not an edit of the first.
- **Solver validation.** DMRG reproduces ED gaps and ground energies to rel ≤ 2e-13 on 4 cases (`mps/validation_vs_ed.json`). The largest ground-state entanglement entropy is 0.83, so chi=120 has a lot of headroom.
- **Bond-dimension check (chi 160 vs 120 at N=192, x=64):** the gap changes by -4.8e-8 (m=0) and -1.8e-8 (m/g=0.25). Truncation is negligible next to the extrapolation spread.

## What this does and does not establish for GRUT

- **It does:** provide a successful spectrum benchmark and a finite interacting realization of the detector-history rule. The frozen spectrum fits pass their ±3% windows. The selected Gaussian description misses finite massive-model histories, including the exploratory matched-variance controls in the Codex review.
- **Matching remains open:** these results do not establish equation (14). The spectrum and history grids are different, the history ladder is not exactly fixed volume, and its Gaussian reference uses the lattice's own mean and covariance. Independent interacting-continuum histories with matched parameters, boundaries and instruments remain necessary. Failure of this Gaussian summary does not exclude every other reduced description.
- **It does not:** This is standard lattice gauge theory placed in GRUT's framework. It does not derive the Schwinger model from GRUT's access or closure principles, and it does not address chirality, 3+1D, or gravity. The open obligations listed in `HANDOFF.md` are unchanged.

## Review and replay scope

The Codex review reran four solver checks, independently checked finite Hamiltonian matrices, recalculated the frozen spectrum fits, and reproduced one production point. It did not rerun the complete production grid or independently implement DMRG. See [the review](../grut_qed_review_20261107/REVIEW.md).

From a clean checkout with NumPy, SciPy and `physics-tenpy==1.1.1` installed, run from any directory, replacing `/path/to/TestingGRUT` with the checkout path:

```bash
run_dir=$(mktemp -d)
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python /path/to/TestingGRUT/research/reality_framework/grut_qed_detector_histories_20261107/mps/validate_vs_ed.py --output "$run_dir/validation.json"
```

The output is separate from archived results. The command returns zero only when the original gate passes; subsets requested with `--case` are explicitly labeled. A supplied `--reference-dir` is recorded with hashes, so a deliberately corrupted-reference check cannot be mistaken for the archived four-case replay. The original `validation_vs_ed.json`, registrations, physics solver and production outputs are unchanged.

# DMRG rerun of P-A and P-B: both pass

The run follows `mps/PREREGISTRATION_AMENDMENT_01.json` (sha256 `0db11998…e18e`), which was sealed in `dd1e107` before any production DMRG. The analysis script was committed in `c4e1ddf` before any results were read. Targets and ±3% windows are the same as in the original `ed/PREREGISTRATION.json`.

| Prediction | Target | DMRG primary | Error | Variant spread (12 fits) | Verdict |
|---|---:|---:|---:|---|---|
| P-A, vector mass at m=0 | 0.5641896 (exact, \(1/\sqrt\pi\)) | 0.5641906 | +1.7e-6 | 0.564179–0.564588 (target inside) | PASS |
| P-B holdout, binding energy at m/g=0.25 | 0.51918 (Byrnes et al. 2002) | 0.5192320 | +1.0e-4 | 0.519207–0.520345 (target 0.005% below lower edge) | PASS |

## How to read it

- **Precision.** The spread is 0.07% at m=0 and 0.2% at m/g=0.25. The agreement of the P-A primary to 2e-6 is partly fortunate and should not be quoted as the method's precision. The meaningful statement is that every variant lies within 0.25% of its target.
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

- **It does:** GRUT's detector-history setup has now been tested on interacting gauge matter with known continuum answers. The same lattice that passes the history tests P-C/P-D/P-E reproduces the bound-state spectrum to well under 1%. The matching condition (14) holds in this test, for 1+1D QED.
- **It does not:** This is standard lattice gauge theory placed in GRUT's framework. It does not derive the Schwinger model from GRUT's access or closure principles, and it does not address chirality, 3+1D, or gravity. The open obligations listed in `HANDOFF.md` are unchanged.

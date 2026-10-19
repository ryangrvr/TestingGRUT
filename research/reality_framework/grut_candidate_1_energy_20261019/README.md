# Candidate 1 quantum energy packet

Start with [REPORT.md](REPORT.md). This packet continues the frozen Candidate 1
scalar/Einstein EFT on the previous fluctuations branch; no historical candidate
or earlier result is overwritten. It computes one scalar-loop energy and a
specific free-carrier comparison, not a lifetime or universal theory.

## Run the verification

Python 3.12, NumPy and SciPy are sufficient. Executed versions are Python
3.12.14, NumPy 2.3.5 and SciPy 1.17.0. In this directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
```

This rechecks stored evidence and performs fresh canonical-toy, continuum-vacuum
and bubble-kernel calculations. It does not pretend to rerun the entire physical
integral. The checked-in audit is a separate full, stricter replay of that integral.

## Recalculate the physical results

From this directory, with the original sibling candidate/fluctuation packets
present in the same repository checkout:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python graphs.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python spectral.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python tail_refinement.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python integrate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
```

These commands replace their own generated output files; run in a disposable
checkout if preserving the published execution bytes is important. They do not
alter the action or input profile. Each large integral has a finite declared
time budget (540 and 600 seconds respectively); a slow machine may require a
clearly labeled budget amendment. Missing points block the binding verdict.
The canonical toy intentionally refuses to overwrite its original output;
`verify.py` reruns it in a temporary directory and compares the result.

Expected integrated values, not exact-bit platform requirements:
`K1/(gv)≈2.79528`, numerical allowance ±.11;
`DeltaE_freeX/(gv)≈-3553.70573` at integer N39944.
The allowance is not an exact-profile, all-orders or empirical uncertainty.

## Evidence map

| Files | Purpose |
|---|---|
| `ENERGY_DERIVATION.md`, `theory_toy.py`, `theory_*json` | Fixed-charge and complete-vacuum normalization |
| `graphs.py`, `graphs_*` | Matched MSbar Born restoration and transform controls |
| `spectral.py`, `spectral_*` | Original continuum determinant and point records |
| `audit.py`, `audit_*` | Full stricter replay, FD method and endpoint controls |
| `tail_refinement.py`, `tail_refinement_*` | Sharper weighted angular bound |
| `integrate.py`, `binding_*` | Integer-charge BVP and like-order pole comparison |
| `SPECTRAL_NOTES.md`, `CHALLENGE.md` | Derivations, actual audit scope and limitations |
| `verify.py`, `VERIFICATION.json` | Executed checks and fresh limited recalculations |
| `NEXT_CALCULATION.md` | Bounded emission/equal-split test, not yet executed |

The manifest hashes this packet's final files and identifies immutable parent
inputs. It is a byte-integrity record, not a physical validation certificate.

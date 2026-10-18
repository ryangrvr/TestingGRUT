# Reproducible leading-regime calculations

Tested environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 on Linux. No additional packages were installed during this construction. Versions are pinned in `requirements.txt`; the first two libraries are sufficient for the physical calculations, while Matplotlib produces the figure.

In a fresh extracted copy, from the package root:

```bash
python code/run_all.py
```

This executes sequentially and terminates. Equivalent component commands:

```bash
python code/qball.py
python code/results_additional_checks.py
python code/reporter.py
python code/verify.py
python code/plot_regime.py
```

Runs overwrite generated outputs in that copy. Preserve the distributed execution records if comparing to the original run. Runtime depends on hardware and eigensolver behavior; the main scalar run took approximately 90 seconds and the original reporter calculation about 4 seconds before the extra packet grids were added. No indefinite execution is intended.

Reproduction means agreement of the numerical observables within the stated convergence limits. Runtime fields, timestamps, adaptive meshes, eigensolver roundoff and compressed-file metadata need not be byte-identical. The final manifest authenticates the distributed snapshot, not every future output file.

## Files and provenance

| Files | Meaning |
|---|---|
| `qball.py` | Fixed-preparation BVP refinements, scalar energy/charge, constrained amplitude spectrum, phase spectrum, stress and gravity |
| `results.json`, `results_prescription.json` | Full final results and physical/numerical prescription; supplemental diagnostics are labeled |
| `results_profile.npz`, `results_profile.csv` | Solved background for reporter/plot calculations; NPZ stores `rho,y=(F,Fprime,S,Sprime),Omega` |
| `results_execution.log`, `results_notes.md` | Main execution record and calculator's scope notes |
| `results_additional_checks.py/json` | Saved-profile recheck, phase refinements and unconstrained radial negative-mode diagnostic |
| `results_preparation_search.py/json`, `results_preparation_seed.npz` | Successful fixed-frequency bag-seed search and supplied preparation |
| `results_initial_*`, `results_gaussian_*` | Six failed preliminary BVP attempts; these are retained failures, not successful states |
| `variational_seed.py`, `variational_seed_results.json` | Exploratory Gaussian trial-function initial-guess planning; not a stationary solution or stability proof |
| `reporter.py`, `reporter_prescription.json` | Three declared monochromatic momenta; two independent phase methods; three packet grids and ideal region/arm instrument |
| `reporter_results.json`, `reporter_execution.log` | Phase/Born-bound and packet outputs, with profile hash and environment |
| `verify.py`, `verification.json` | 25 acceptance checks on equations, identities, convergence, spectra and probabilities; not an empirical or global-stability test |
| `plot_regime.py` | Scientific profile/stress/metric/phase plot from saved outputs |
| `run_all.py` | Sequential finite reproduction wrapper |

The scalar solve fixes `F(0)` and solves for `Omega`; charge is an output. The subsequent Hessian tests variations at that computed fixed charge and correctly include the positive rank-one radial term. These are different tasks, not an inconsistency in the preparation specification.

The main script reuses the saved numerical preparation. If that seed is absent, it executes the included preparation search, retaining its search outputs. Earlier failed-start scripts are not re-executed by `run_all.py`; their configurations, outputs and logs remain in the archive. The Gaussian variational script can be run independently, but its trial-function output must never be substituted for the BVP solution.

## What reproduction does and does not reproduce

It reproduces a stationary classical scalar profile, local fixed-charge energetic spectra, full-stress first-order Einstein response, classical linear reporter scattering, and a supplied positive-frequency one-particle/ideal arm readout. The finite packet uses a radial L2 representation and unitary matrix evolution; norm conservation alone does not establish continuum accuracy. Its grid sweep is printed explicitly.

It does not simulate a physical source/detector, all Standard Model dynamics, a nonlinear gravitating lump, arbitrary matter formation, a complete quantum target/probe system, quantum lifetime, every fragmentation channel or experimental data. All physical conclusions retain the limits stated in the theory documents.

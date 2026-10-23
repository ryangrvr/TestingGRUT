# Reproduce the completed inexpensive screen

Run from this TestingGRUT commit. The runtime used Python 3.12, NumPy 2.3.5,
SciPy 1.17.0 and Matplotlib 3.10.8. No system packages were installed.
Only the pinned hipsofcobra source is executed in temporary scratch.
Its generated tables are included; the external public toy archive is
approximately 33 MB compressed and needs substantially more working memory.

```sh
TASK_PACKET=research/reality_framework/grut_portal_debris_window_20261023
TASK_CACHE=/tmp/grut-portal-window-inputs
export OPENBLAS_NUM_THREADS=1
python "$TASK_PACKET/code/make_manifest.py" --check
python "$TASK_PACKET/code/fetch_inputs.py" "$TASK_CACHE" --regenerate-hadronic
python "$TASK_PACKET/code/profiles.py"
python "$TASK_PACKET/code/radial.py"
python "$TASK_PACKET/code/scan.py"
python "$TASK_PACKET/code/projection.py" "$TASK_CACHE"
python "$TASK_PACKET/code/heat.py"
python "$TASK_PACKET/code/fixed_yield.py" "$TASK_CACHE"
python "$TASK_PACKET/code/verify.py" --cache "$TASK_CACHE"
python "$TASK_PACKET/code/plot.py"
```

`profiles.py` attempts 17 finite-Q profiles and retains the old-charge
continuation failure. `radial.py` adds its original old-charge refinement
and six finite-grid radial amplitude Hessians. `scan.py` writes 1021
discrete rows and eight selected points. `projection.py` uses all 1000
background-only public toys for those eight, not a discovery forecast.
`heat.py` reuses the parent's conditional redshift/heat kernel.

`fixed_yield.py` adds the Q=200000 profile with larger-domain refinement,
then calibrates v for two test populations at the original supplied YQ.
Its charges are not predicted. `verify.py` separately reconstructs raw
profile integrals, invariant two-body phase-space and population budgets,
checks detector quadrature, and refines the public-toy strength grid.
It does not run a new BBN, inelastic lifetime, phase-history or LZ fit.

The six numerical programs were also run in a fresh isolated packet copy;
1212 numerical leaves agreed with the recorded outputs at stated
floating-point comparison tolerances. `REPLAY_RESULTS.json` records the
counts. It is same-author reproduction, not independent review. Source
hashes and regenerated hadronic tables are separate checks.
Run a numerical replay in a copy of the delivered packet if preserving its
exact hashes matters: floating-point outputs can differ slightly across
platforms even when scientific comparisons agree. The initial manifest
check verifies the delivered bytes before any output is regenerated.

After an intentional packet change, regenerate the manifest only after
the numerical outputs and documentation are final:

```sh
python "$TASK_PACKET/code/make_manifest.py"
python "$TASK_PACKET/code/make_manifest.py" --check
```

The manifest includes consumed parent inputs. Larger old archive results
retain their published manifests and scope rather than being recertified
here. A successful replay establishes reproducibility of the stated
calculation, not existence of the proposed fields or cosmological viability.

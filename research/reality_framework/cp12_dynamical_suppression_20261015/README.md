# CP-12: dynamical suppression under the full frozen law

Continuation of `768c46c96eae404a527d9cfee8f1b09abeb61a13`, on the isolated
branch `codex/grut-dynamical-suppression-20261015`. Canonical GRUT-RAI,
Stage3 and PR12 are not changed.

Read [REPORT.md](REPORT.md) for numerical findings,
[DERIVATION.md](DERIVATION.md) for definitions and error control, and
[HANDOFF.md](HANDOFF.md) for the next builder/reviewer decision.

The packet tests Pauli availability, diagonal detuning and coherent
interference while evolving the **full** original Hamiltonian. It also
constructs larger six-site numerical bases with every omitted physical
transition retained in the error residual. No hopping term is removed in
these new evolutions. The algebraic supplied-state controls are explicitly
separate from the neutral-origin formation experiments.

## Reproduce

Run from the repository root with Python 3, NumPy, SciPy, Matplotlib and g++.
The recorded environment used Python 3.12, NumPy 2.3.5 and SciPy 1.17.0.
The inherited symmetry helper compiles locally. Large caches are optional
and are not research inputs; the scripts can regenerate them from the
unmodified physical destination routine.

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GRUT_SHORTCUT_CACHE=/tmp/grut-shortcut-cache
export GRUT_SUPPRESSION_CACHE=/tmp/grut-suppression-cache
packet=research/reality_framework/cp12_dynamical_suppression_20261015

python "$packet/checks.py"
python "$packet/pauli_algebra.py"
python "$packet/neutral_exterior.py"
python "$packet/run.py" --n 5 --budget 50 --waits .5,1,2 --lags .125,.25,.5 --out "$packet/N5_B50.json"
python "$packet/rare_pauli.py"
python "$packet/core_response.py"
```

Build the six-site extension in this order. The priority is only a basis
allocation proxy; its coverage is not a probability or an error bound.

```bash
python "$packet/priority.py" --out "$GRUT_SUPPRESSION_CACHE/n6_priority.npz"
python "$packet/boundary.py" --priority "$GRUT_SUPPRESSION_CACHE/n6_priority.npz" --count 225000 --name N6_ADAPT225K
python "$packet/run.py" --n 6 --budget 12 --model N6_ADAPT225K --waits .5,1,2 --pilot-only --out "$packet/N6_ADAPT225K_PILOT.json"
python "$packet/boundary.py" --parent-model N6_ADAPT225K --first 480212 --name N6_SHELL1
python "$packet/boundary.py" --parent-model N6_SHELL1 --all --name N6_SHELL2
python "$packet/run.py" --n 6 --budget 12 --model N6_SHELL2 --waits .5,1 --lags .125,.25,.5 --out "$packet/N6_SHELL2.json"
python "$packet/boundary_checks.py"
python "$packet/plot.py"
python "$packet/verify.py"
```

Run large matrix construction and propagation **serially**. An unfinished
intermediate propagation concurrent with construction, and an initial
explicit-generator run, exceeded this workspace's roughly 10 GiB memory
allowance. Their incomplete trajectories are not used as numerical
evidence. The final six-site implementation reuses real $H$
through a LinearOperator and releases basis labels after constructing all
observables. An independent calculation agrees with the explicit generator
to about $10^{-15}$.
The original complete budget-12 closure, priority extension and successive
shells remain reproducible numerical approximations to the same full law.
The uncertainty intervals, rather than visually stable central estimates,
determine what can be concluded.

## Implementation map

| File | Purpose |
|---|---|
| `run.py` | Original-neutral preparation, common quarter-threshold event, full-H action integrals and bounded mechanism diagnostics |
| `mechanisms.py` | Physical-basis interference references, neutral endpoint availability, capped diagonal detuning; excludes outgoing rewiring from hopping action |
| `action_integrals.py` | Squared-action midpoint quadrature, analytic remainders, outgoing-residual state bounds and full-law action intervals |
| `hopping_bound.py` | Exact integer certificate for global filling-sector bounds on both hopping actions |
| `linear_propagator.py` | Real-H LinearOperator and exact adjoint, avoiding duplicate large sparse generators |
| `resources.py` | Complete binary-record outcome entropy derived from the formation-probability interval |
| `priority.py` | Deterministic outgoing-fluence proxy used only for numerical basis allocation |
| `boundary.py` | Exact signed-orbit extensions and every omitted outgoing transition |
| `boundary_checks.py` | Sampled large-basis columns independently regenerated with the original Python canonicalizer |
| `checks.py` | Independent labelled physical comparison, spectral quadrature check and complete small adaptive-column check |
| `pauli_algebra.py` | Exact rational blocking/unlocking and extensive-denominator controls; supplied states only |
| `neutral_exterior.py` | Coherent physical-basis check of the exact initial neutral-exterior action identity |
| `rare_pauli.py` | Exploratory neutral-origin, rare filled-ring event; formation cost and exact conditional shape |
| `core_response.py` | Auxiliary N4,Q5 full-law response at fixed recorded-ring distances; complete record resources |
| `plot.py` | Research figure from saved results |
| `verify.py` | Numerical assertions, frozen-source hashes and content manifest |

All numerical `PASS` fields mean the executable checks completed, not that
CP-12 passed the physical locality criterion. Ordinary floating-point
intervals include analytic truncation and quadrature bounds plus a reporting
allowance; they are not outward-rounded formal certificates.

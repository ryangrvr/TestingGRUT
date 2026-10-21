# Reproduction and provenance

Start from this packet's Git commit. Public pre-execution freezes are:

- `188b92b9375df468802cb064c8a288bcc38f956f`: original numerical protocol.
- `254254ea8d35c5ea0355c28b4fce601243d56bf5`: estimator/initialization clarification.
- `8ed29f00b402054a5d88cbb796e1dfb1cb8cdf2b`: reduced timestep after failed orbit control.
- `998ddd8b1252e3ecd7cab8be4dc99044574f5faa`: validated exact source hashes and main run list.
- `a9688dc825484eeee61ea11d930880734ea55313`: exploratory Gaussian protocol after early sparse results.
- `ed998463ce01aa1f4d352e04be59cdd9ae1f79a1`: validated Gaussian wrapper hashes before Gaussian production.

Source read: Claude `c7b8480a49e3f0b8737e800396a30f33f6337981`, report and issue #13
comment 6314157022; charged-exit/bias notes and selected code and verification records.
The complete twelve-module calculation campaign was not rerun. The earlier Codex
Cartesian solver was read and the same normalized action/stencil were independently
extended here. No Claude Python integrator was imported or copied; the conformal-variable
mathematical method is shared. The Gaussian wrapper reuses this packet's frozen C++
integrator, and is not a second independent integrator.

To inspect the recorded diagnostics without rerunning production:

```bash
cd research/reality_framework/grut_expanding_charge_20261021
python archive.py unpack
OPENBLAS_NUM_THREADS=1 python verify.py
```

Dependencies: g++ with C++17/OpenMP, Python, NumPy and SciPy. Matplotlib is used only
for the figure. Exact session versions, flags and hashes are in ENVIRONMENT.json,
SOURCE_FREEZE.json, GAUSSIAN_SOURCE_FREEZE.json and each raw run configuration.
There is no GPU requirement and no `-march=native` or fast-math build.

Controls and frozen suite:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=8 python run.py validate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=8 python run.py production
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=8 python gaussian_run.py validate
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=8 python gaussian_run.py production
OPENBLAS_NUM_THREADS=1 python analyze.py
```

The launch caps depend on machine speed, so replaying the suite on a faster machine
can launch jobs that were not executed in this session. RUN_STATUS.json and
GAUSSIAN_STATUS.json explicitly distinguish those jobs. For the exact executed list,
use the individual argument dictionaries in the archived configuration files or
`python reproduce.py --list`, then select a run with `--run NAME`.

The original run's commands contain this session's absolute paths for provenance.
The drivers derive current paths from their own location, so replay does not require
those paths. Runtime and metadata hashes need not match on replay. Parallel floating
point reductions can change last digits of integrals; field initialization is
deterministic. Compare quantities and stated convergence, not just file hashes.

No full evolving 3D field dumps are archived. All scalar histories, operational
topology histories, mode tables, configurations and run metadata are losslessly
archived. Central f32 slices were visualization aids; their hashes are recorded
but their arrays are not included in the small archive. Replaying a run recreates
its slice. Plot.py needs those slices to regenerate the lower image panels.

The two failed pre-production records are retained. The first had a nondeterministic
normalization mismatch and failed the coarse orbit limit; the second still failed
the coarse orbit limit. The corrected source and reduced timestep passed before
production. These are internal numerical controls, not independent outside scientific
replication or empirical observations.

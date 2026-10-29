# Reproduce the kinetic packet

Use Python with NumPy, SciPy and Matplotlib, plus a C++17 compiler. Exact
versions are in ENVIRONMENT.json. Set BLAS threads to one for reproducible
resource use. Run from the TestingGRUT root at this packet's commit.

```bash
export OPENBLAS_NUM_THREADS=1
export GRUT_M1_KINETICS_WORK=/tmp/grut_m1_kinetics_reproduction
python research/reality_framework/grut_m1_kinetics_20261028/code/source.py
python research/reality_framework/grut_m1_kinetics_20261028/code/core.py
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode core --grid translated --shift-bins 64 --emax 8 --tag translated_core_64
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode resonant --grid translated --shift-bins 64 --emax 8 --tag translated_resonant_64
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode full --grid translated --shift-bins 32 --emax 8 --angles 12 --tag translated_full_32
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode full --grid translated --shift-bins 64 --emax 8 --angles 12 --tag translated_full_64
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode full --grid translated --shift-bins 96 --emax 8 --angles 12 --tag translated_full_96
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode full --grid translated --shift-bins 64 --emax 12 --angles 12 --tag translated_full_64_range
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode full --grid translated --shift-bins 64 --emax 8 --angles 16 --tag translated_full_64_angles
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode full --grid aligned --shift-bins 64 --emax 8 --angles 12 --tag aligned_full_64
python research/reality_framework/grut_m1_kinetics_20261028/code/transport.py --mode full --grid aligned --shift-bins 64 --emax 8 --angles 12 --rtol 2e-9 --tag aligned_full_64_tight
python research/reality_framework/grut_m1_kinetics_20261028/code/physics_checks.py
python research/reality_framework/grut_m1_kinetics_20261028/code/verify.py
python research/reality_framework/grut_m1_kinetics_20261028/code/figure.py
python research/reality_framework/grut_m1_kinetics_20261028/code/audit.py
```

Use an empty work directory for fresh source and core-kernel reproduction.
Core-kernel cache keys include the profile and implementing source hashes.
The inherited packet files and engine are read, not overwritten.

Each transport JSON records its full configuration. Its matching NPZ contains
the energy grid, optical times, all species populations, the omission counter,
and accepted ODE step times. Earlier unaligned, uniform and exploratory long
runs are controls, not the primary numerical conclusion. Their exact command
options are retained in their JSON configurations. The initial endpoint-bug
outputs are preserved under attempts/floor_rejection/ and cannot be regenerated
with the corrected source without reverting that historical algorithm.

physics_checks.py additionally reads the retained long-run control if present;
that row is optional and is not used to pass a primary-history test. The report
designates the translated graded 96-node-per-interval result as primary.

The complete source+core+refinement sequence took minutes in the recorded
environment. No cosmological simulation or external AI service is required to
run it. Agreement with these numbers is method replication, not observational
confirmation of GRUT-M1.

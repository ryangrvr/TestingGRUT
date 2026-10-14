# CP-12 relational record gate

Read [REPORT.md](REPORT.md) for the result and [DERIVATION.md](DERIVATION.md) for the definitions and bounds.

This packet analyzes the existing states from `8ece3e2` after the exchangeability audit `a351fea`. It adds no new physical Hamiltonian, preparation, vertex size, source, recording operation or longer evolution time. Previously saved observable files did not contain wavefunctions; `analyze.py` reconstructs the same projected states and checks the replay against every saved count. Optional snapshots avoid further replay during analysis. Private cache pickles are only a speed optimization and are not required.

Install the pinned dependencies with Python 3.11 or later. From this directory, a minimal reproduction is:

```bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python analyze.py --n 4 --q 4 --budget 256 --snapshot-dir /tmp/cp12-relational-waves
python checks.py
python energy_control.py
python verify.py
```

The complete already-tested cases are:

```bash
python analyze.py --n 4 --q 4 --budget 128 --times 0 1 3 5 10
python analyze.py --n 5 --q 5 --budget 128 --snapshot-dir /tmp/cp12-relational-waves
python analyze.py --n 6 --q 6 --budget 128 --snapshot-dir /tmp/cp12-relational-waves
python analyze.py --n 6 --q 12 --budget 10 --times 0 0.4 0.6 1 --snapshot-dir /tmp/cp12-relational-waves
python refine.py --snapshot-dir /tmp/cp12-relational-waves
python mechanism.py --snapshot-dir /tmp/cp12-relational-waves
python shortcut_controls.py --snapshot-dir /tmp/cp12-relational-waves
python plot.py
python verify.py
```

Set `CP12_FORMATION_CACHE` to an optional trusted private cache directory to reuse matrices from the previous packet. Otherwise the original model builder reconstructs them. Never load untrusted pickles. `refine.py` uses auxiliary backward vectors to estimate errors; this is not a physical reversal or intervention.

`verify.py` freshly reruns the independent small-graph and physical-orbit checks, the N4 analysis, and the mean-energy control in temporary output directories. It validates the larger saved calculations and their hashes; it does not silently rerun the larger replay or dense backward refinement. `verify.py --record-manifest` intentionally replaces the hash manifest after a reviewed source/output update. No external files are downloaded by the verifier.

The Gauss/resource references have a restricted support and an explicit identity extension outside it. They match diagonal resources, not full Hamiltonian energy. The separate physical N4 spectral and maximum-entropy references match mean energy but are mixed and do not match each instantaneous resource distribution. These different controls must not be conflated.

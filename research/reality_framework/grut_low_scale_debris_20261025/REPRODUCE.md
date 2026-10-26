# Reproduction

The main production needs Python 3.12, NumPy, SciPy, Matplotlib and GCC.
Temporary release output is outside the repository, with a configurable
GRUT_LOW_DEBRIS_WORKDIR (default /tmp/grut_low_scale_debris_20261025).

From this packet's code directory:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python build.py
python releases.py
python decays.py
python vector.py
python benchmark.py
python finish_ledger.py
python report.py
python controls.py
python manifest.py
```

The source release engine is inherited explicitly from the parent packet.
Four new full releases check long windows and tighter tolerances; all 96
parent ensemble endpoints are independently projected. Main relic results
come from the newly written vectorized Bose integrals and a positive
implicit exponential midpoint update, not imports of reviewer production
code. Steps 260/520/1040, quadratures 96/128, the late endpoint and a
tenfold reaction-rate range are recorded. Initial numerical issues with
near-zero charge subtraction and endpoint brackets were corrected using
algebraically stable variables; no failing physical trajectory was dropped.

For the optional **one-case reviewer-code replay**, first extract a
read-only archive of commit f1aa7df968cb3dc90ca32f4f32aeb0a819ddfce4 outside
the working tree, then run:

```bash
python reference_replay.py /path/to/review_archive_root
```

The result is saved in inputs/REVIEW_BENCHMARK_REPLAY.json. It uses the
reviewer's original adaptive thermodynamic integrals, 260-step endpoint
update and coarse entropy table. It is a cross-check of this benchmark,
not outside replication. Run it before controls.py to include that check.

The main calculation replaces the review entropy table's nonmonotone
g_s switch with a separate continuous entropy-degree interpolation; the
legacy replay is retained to disclose its numerical effect. A source-only
case is also continued below the old 1-keV endpoint. All physical statuses
remain separate from numerical checks.

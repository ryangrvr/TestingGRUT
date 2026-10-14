# CP-12 flux-loop gate

Read [REPORT.md](REPORT.md) and [DERIVATION.md](DERIVATION.md). This additive packet uses the unchanged frozen law, the supplied neutral ray, existing evolved states and the physical N4 energy references from `9d0f2ef`. The parent is `fe90e7f`.

The primary populations are unmeasured-state diagnostics. Persistence and density-response files clearly distinguish a hypothetical initial loop filter, unforced two-time correlators, source kicks, and final circulation conditioning. These diagnostic operations are not an autonomous topology selector.

Install the pinned requirements with Python 3.11 or later. From this directory:

```bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python onset_exact.py
python controls.py
python references.py
python early_time.py
python verify.py
```

`verify.py` runs fresh independent small controls, reconstructs the N4 snapshots and checks the main same-state diagnostic. It validates larger saved outputs and hashes without silently repeating the N5/N6 replay, dense dual estimate or all mixed-state persistence curves.

To reproduce all population files, first obtain the snapshots using the preceding packet's analyzer. Pickles are an optional trusted private speed cache; none is required for reproduction.

```bash
python ../cp12_relational_record_gate/analyze.py --n 4 --q 4 --budget 256 --snapshot-dir /tmp/cp12-loop-waves
python ../cp12_relational_record_gate/analyze.py --n 5 --q 5 --budget 128 --snapshot-dir /tmp/cp12-loop-waves
python ../cp12_relational_record_gate/analyze.py --n 6 --q 6 --budget 128 --snapshot-dir /tmp/cp12-loop-waves
python ../cp12_relational_record_gate/analyze.py --n 6 --q 12 --budget 10 --times 0 0.4 0.6 1 --snapshot-dir /tmp/cp12-loop-waves
python populations.py --snapshot-dir /tmp/cp12-loop-waves
python onset_uniform.py
python persistence.py --snapshot-dir /tmp/cp12-loop-waves
python circulation_survival.py --snapshot-dir /tmp/cp12-loop-waves
python measurement_controls.py --snapshot-dir /tmp/cp12-loop-waves
python transport.py --snapshot-dir /tmp/cp12-loop-waves
python refine_dense.py --snapshot-dir /tmp/cp12-loop-waves
python plot.py
python verify.py
```

Setting `CP12_FORMATION_CACHE` reuses matrices from earlier work. Do not load untrusted pickles. The `source=True` builder flag used internally selects a computational stabilizer basis; it inserts no physical source or Hamiltonian term. Explicit source kicks occur only in `transport.py`.

Electric circulation is checked on the full incidence graph. It can include hole edges; it must not be equated with occupied-edge support. The Q=N spanning circulating networks are simple rings, not a three-dimensional scaling family. Projected zeros in the dense budget10 calculation are accompanied by full-law upper bounds.

`verify.py --record-manifest` intentionally replaces the manifest after reviewed source/output changes. The numerical bounds do not include a formal interval-arithmetic roundoff certificate. Exact rational onset coefficients are separately identified.

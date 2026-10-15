# CP-12 shortcut dynamics — builder continuation

Read `REPORT.md` for the finding, `DERIVATION.md` for the operators and bounds,
and `HANDOFF.md` for the reviewer continuation. Parent: `aedf492`.

Python 3.12, NumPy, SciPy and a C++17 compiler are sufficient. The inherited
integer canonicalizer compiles its helper automatically. No result depends on
a local cache or an unavailable trajectory file. Optional basis caches can be
enabled by setting `GRUT_SHORTCUT_CACHE` to a writable directory.

From the repository root, reproduce the evidence with:

```bash
python research/reality_framework/shortcut_locality_builder_20261015/shortcut_control.py > research/reality_framework/cp12_shortcut_dynamics_20261015/PARENT_TOY_REPRODUCTION.json
python research/reality_framework/cp12_shortcut_dynamics_20261015/checks.py
python research/reality_framework/cp12_shortcut_dynamics_20261015/transfer_algebra.py
python research/reality_framework/cp12_shortcut_dynamics_20261015/run_split.py --n 5 --budget 50 --times .25,.5,1,2,3 --wait 2 --lags .125,.25,.5 --out research/reality_framework/cp12_shortcut_dynamics_20261015/SPLIT_N5_B50.json
python research/reality_framework/cp12_shortcut_dynamics_20261015/run_split.py --n 6 --budget 12 --times .25,.5,.75,1 --wait .5 --lags .125,.25 --out research/reality_framework/cp12_shortcut_dynamics_20261015/SPLIT_N6_B12.json
python research/reality_framework/cp12_shortcut_dynamics_20261015/matter_response.py --n 5 --budget 50 --wait 2 --lags .125,.25,.5 --out research/reality_framework/cp12_shortcut_dynamics_20261015/RESPONSE_N5_B50.json
python research/reality_framework/cp12_shortcut_dynamics_20261015/recorded_transport.py --budget 128 --wait 5 --out research/reality_framework/cp12_shortcut_dynamics_20261015/TRANSPORT_N4_Q5_B128.json
python research/reality_framework/cp12_shortcut_dynamics_20261015/recorded_transport.py --budget 256 --wait 5 --out research/reality_framework/cp12_shortcut_dynamics_20261015/TRANSPORT_N4_Q5_B256.json
python research/reality_framework/cp12_shortcut_dynamics_20261015/plot.py
python research/reality_framework/cp12_shortcut_dynamics_20261015/verify.py --write-manifest
python research/reality_framework/cp12_shortcut_dynamics_20261015/verify.py
```

`OPENBLAS_NUM_THREADS=1` was used for predictable memory/CPU usage in this
execution. Fresh N6 closure and the N5 source-stabilizer closure take minutes.
`verify.py` checks saved evidence and integrity; it does not rerun the complete
propagations. The complete commands above are the numerical reproduction.
Timing fields vary by machine; `--write-manifest` records the independently
checked reproduction's new file hashes. Use plain `verify.py` to verify the
unchanged published packet without replacing its manifest.

`operators.py` constructs both the exact signed permutation quotient and the
independent labelled physical basis. `integrals.py` evaluates multiple
outgoing-action integrals without retaining a full time grid. `checks.py`
compares these with the inherited independent implementations.

All full-law results preserve the frozen electric, mass, hopping and rewiring
parameters. Only explicitly named comparison evolutions remove hopping.
Neither the four-site auxiliary sector nor the five-/six-site dense hosts
establish a large-distance locality limit.

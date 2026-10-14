# Distance before and after a source intervention

Use a complete checkout containing the parent CP-12 packets. Python 3, NumPy and SciPy are required. All new records and code are additive; no frozen physics source is modified.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/cp12_distance_history/verify.py
```

The verifier runs the parent checks, fresh ordered-pair/full-physical-matrix controls, exact edge-pool incidence checks, matched-state attribution controls, and finite-pulse extraction. It validates saved evolution and state-exposure records separately; it does not silently repeat the six expensive evolutions or the state-exposure audit.

To regenerate the main recorded/unrecorded comparisons:

```bash
P=research/reality_framework/cp12_distance_history
OPENBLAS_NUM_THREADS=1 python "$P/experiment.py" --n 4 --budget 10 --lag .1
OPENBLAS_NUM_THREADS=1 python "$P/experiment.py" --n 4 --budget 10 --lag .3
OPENBLAS_NUM_THREADS=1 python "$P/experiment.py" --n 5 --budget 10 --lag .1
OPENBLAS_NUM_THREADS=1 python "$P/experiment.py" --n 5 --budget 10 --lag .3
OPENBLAS_NUM_THREADS=1 python "$P/experiment.py" --n 6 --budget 8 --lag .1
OPENBLAS_NUM_THREADS=1 python "$P/experiment.py" --n 6 --budget 8 --lag .3
OPENBLAS_NUM_THREADS=1 python "$P/attribution.py"
OPENBLAS_NUM_THREADS=1 python "$P/state_coupling.py"
OPENBLAS_NUM_THREADS=1 python "$P/pulse.py"
OPENBLAS_NUM_THREADS=1 python "$P/edge_speed.py"
python "$P/summarize.py"
```

Q defaults to N−1 and the waiting time to 0.1. These sparse sectors are auxiliary propagation tests; connected members are trees. They do not replace the higher-density gauge-capacity requirement. The pool theorem is independently stated for general Q.

The six-event source–receiver quotient at budget eight has 53,169 vectors and roughly one million Hamiltonian entries. Recorded runs took several minutes on the host; runtime and RAM vary. Optional CP12_FORMATION_CACHE points to a TRUSTED local scratch directory holding disposable pickle caches. Omit it for uncached reproduction. Do not load third-party pickle files.

The main JSONs contain the pre-distance class weights, conditional response intervals, joint signed histories, positive distance-transition weights, recorded and unrecorded post-distance responses, and recording energy/entropy. The disconnected class is N numerically. Labels such as class=6 at N=6 mean disconnected, not graph distance six.

The six-event lag-0.1 record predates an optional saved norm check for the coherent sum of branches; its history sum rules are still verified. Later records include that additional field. Grouped propagation/refactoring can change final floating-point digits and runtimes; the operator definitions and analytic bounds are unchanged. The small complete-matrix controls test both the ordered quotient and purification/instrument rules.

All error bounds are analytic bounds on omitted-state effects and residual quadrature, evaluated numerically. They are not formal interval-arithmetic bounds on floating-point roundoff. An independent finite-pulse computation checks the five-event signal to much smaller numerical differences than the reported response. Tiny values with equally tiny truncation bounds should not be presented as machine-certified intervals.

REPORT.md states the locality ruling and limits. DERIVATION.md defines the measurement instrument, proves the bounds, separates recording backaction from distance changes, and derives the exact collective-pool coupling norm. --past-only is a computational option for exploratory runs and does not generate the full saved comparison schema.

# Reproduction

Use Python with NumPy and SciPy. The reference run used Python 3.12.14,
NumPy 2.3.5 and SciPy 1.17.0. Run these commands in this directory.
The timeout is per solver, and thread counts are fixed for the run.

~~~bash
replay_dir=$(mktemp -d)
timeout 120s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python i3_calculate.py --output "$replay_dir/production.json"
timeout 120s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python critical_replay.py --output "$replay_dir/critic.json"
python compare_replays.py "$replay_dir/production.json" "$replay_dir/critic.json" \
  --output "$replay_dir/crosscheck.json"
python compare_replays.py RESULTS.json "$replay_dir/critic.json" \
  --output "$replay_dir/archive-check.json"
~~~

The production program exponentiates separate labelled sensor and response
blocks, applies the two-state storage kernel, and retains the full
device endpoint distribution for reset. The critical implementation
uniformizes the full 128-state chain from a separately coded transition
enumeration and uses a Poisson tail test. They share only the registered
inputs. The comparison covers 32 rows, all five arms per row, 16 arm
statistics, and two conditional yields per row. Check convergence against
the frozen $10^{-10}$ absolute threshold. The replay loads no external
data and uses no random seed. Results are floating point and can vary
slightly by platform. The study uses analytic probability statements for
optimized policy bounds and binomial tails for its independent scaling
assay.

The reference run agrees within $6.28\times10^{-13}$ in the full
comparison and within $1.36\times10^{-14}$ for bond probability.
The independent omitted mass bound is below $10^{-14}$ per propagation.
Both implementations completed within the individual 120-second caps.
Common clock and apparatus work are admission conditions, distinct
from the computed ideal thermodynamic charges.
[Frozen input](PARAMETERS.json) · [Results](RESULTS.json)
· [Scalar comparison](CROSSCHECK.json)

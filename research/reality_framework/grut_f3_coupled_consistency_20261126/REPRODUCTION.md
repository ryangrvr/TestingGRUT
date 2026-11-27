# Reproduction

The exact numerical checks require Python 3.12 and the
standard library. Run inside this packet directory.
Each implementation has 30 seconds in the frozen budget.

~~~bash
f3_replay_dir=$(mktemp -d)
timeout 30s python replay_f3.py --output "$f3_replay_dir/primary.json"
timeout 30s python critical_replay_f3.py --output "$f3_replay_dir/critical.json"
python compare_f3.py "$f3_replay_dir/primary.json" "$f3_replay_dir/critical.json" \
  --output "$f3_replay_dir/comparison.json"
cmp RESULTS.json "$f3_replay_dir/primary.json"
cmp CRITICAL_RESULTS.json "$f3_replay_dir/critical.json"
cmp CROSSCHECK.json "$f3_replay_dir/comparison.json"
~~~

The primary script evaluates the declared background fractions,
derivative and curvature bounds, symmetric Hessian rank
test, and finite matching certificates. The independent script
uses log derivatives and a separate symmetric-coordinate
contraction and obtains an interval for the exact matter
dispersion. The comparison accounts for the different
coordinate orders for the off-diagonal kinetic entries.
It covers 107 exact scalar values and checks that the
independent interval lies within the primary dispersion bound.
These checks need no random seed, grid or convergence tolerance.

The constraint closure, physical reduction, existence of
the canonical comparison history and the open-branch theorem
are analytical results. The scripts reproduce the frozen
algebraic validation and normalized reality test. They do
not evolve a cosmological trajectory or estimate uncontrolled
quantum corrections. [Primary](RESULTS.json)
· [Critical](CRITICAL_RESULTS.json) · [Cross-comparison](CROSSCHECK.json)

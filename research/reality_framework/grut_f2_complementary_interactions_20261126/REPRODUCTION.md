# Reproduction

The executable checks use Python 3.12 and its standard library.
They use exact rational arithmetic, with a frozen limit of
30 seconds per implementation and 17 sites. Run from this
directory in the repository. The parent F1 directory is required
at the relative path in PARAMETERS.json, with the original
WITNESSES.json content and hash.

~~~bash
f2_replay_dir=$(mktemp -d)
timeout 30s python replay_f2.py --output "$f2_replay_dir/primary.json"
timeout 30s python critical_replay_f2.py --output "$f2_replay_dir/critical.json"
python compare_f2.py "$f2_replay_dir/primary.json" "$f2_replay_dir/critical.json" \
  --output "$f2_replay_dir/comparison.json"
cmp RESULTS.json "$f2_replay_dir/primary.json"
cmp CRITICAL_RESULTS.json "$f2_replay_dir/critical.json"
cmp CROSSCHECK.json "$f2_replay_dir/comparison.json"
~~~

The primary program uses rational Gauss–Jordan inversion and
a physical root of the pair quadratic to reconstruct the springs.
The separate internal implementation uses rational
$LDL^{\mathsf T}$ solves, and independently computes the
certificate constants. The final comparison covers 275 exact
common scalar values, and verifies separately that the distinct
time-tail and joint-response ratio bounds are compatible.

The continuum result is analytic. Its executables replay
time support, profile derivative and volume inequalities,
and constants in the controlled estimates. No grid, trajectory
fit or convergence tolerance is needed. Support equivalences
and the finite-pair inverse theorem are proved in the
accompanying derivations.
[Exact responses](RESULTS.json) · [Critical outputs](CRITICAL_RESULTS.json)
· [Cross-comparison](CROSSCHECK.json)

# Reproduce

From this packet directory, using Python, numpy, scipy and matplotlib:

```sh
OPENBLAS_NUM_THREADS=1 python code/calculate.py
python code/plot.py
python code/verify.py
OPENBLAS_NUM_THREADS=1 python code/verify.py --replay
```

The parent seed must be present at the relative location recorded in
INPUTS.json and have its exact SHA-256. It is part of the same git tree;
do not replace it with a different profile. The calculator imports no old
solver, reads no observations and changes no earlier packet. It solves
five actual integer-charge backgrounds, refines the three main ones,
computes scalar scattering and radial diagnostics, and runs 51 internal
conditions. The exact finite-charge strengthening was after-result work,
recorded in AMENDMENT_01.md, with physical inputs unchanged.

The default verifier is read only. It checks hashes and source identity,
reconstructs energy/charge and the sum rule from the three raw profiles,
uses an independent Numerov finite-difference scattering method, assembles
the finite-mode density matrices explicitly, and independently sums the
binomial decision probabilities. None of these is outside replication.

`--replay` creates a temporary source tree, copies the code, input and hashed
seed, reruns the new calculation, requires complete identical scientific
JSON key sets and compares all eight CSV arrays. Solver evaluation counts
are excluded from physical-leaf comparisons. Environment details are
recorded separately. It deletes the temporary tree afterwards.

For a new draft, render the figure and finish documentation before running
`python code/make_manifest.py`. The manifest includes all packet files
except itself and Python caches. The verifier does not regenerate hashes
or repair altered results. Reproduce published science in a temporary copy
or with `--replay` rather than overwriting the archived output.

No new dependencies were installed. The calculation is local, finite and
exploratory; there is no ongoing background workflow. No old determinant,
lattice, wall history or cosmological likelihood is rerun by these commands.

# Reproduce

From this directory, with Python, numpy, scipy and matplotlib:

```sh
OPENBLAS_NUM_THREADS=1 python code/calculate.py
python code/plot.py
python code/verify.py
python code/verify.py --replay
```

`calculate.py` uses only `INPUTS.json`, writes the scientific results,
environment and CSV files, and runs 42 internal checks. It imports no
previous packet's calculation. The receiver-only Gamma_vis ratios are
declared inputs. A normal calculation does not fetch data or modify any
old scientific packet. `plot.py` renders the static overview.

The default verifier is read only. It checks the manifest, raw mode
numbers/energy, raw Friedmann constraints, a separate retarded radiation
integral and the scientific result dispositions. `--replay` also copies
the source and inputs to a temporary directory, executes a fresh run,
checks every scientific JSON leaf and raw CSV shape/value, then removes
that temporary calculation. It is reproducibility of the same code,
not independent scientific replication.

Solver evaluation counts and environment runtime are computational
diagnostics, not physical leaf comparisons. Complete scientific key
sets and array shapes are required. Unexplained missing output fields
or unrecorded physical input changes are failures.

`MANIFEST.json` is generated last and excludes itself. Regenerating
plots or environment metadata can change hashes. Compare a reproduced
run in a separate checkout or temporary directory; do not replace
published raw files merely to make their hashes agree. The verifier
does not update a manifest to hide a changed output.

No new dependencies were installed for this continuation. The original
three scientific packets, the q review and Candidate 1 are not changed.

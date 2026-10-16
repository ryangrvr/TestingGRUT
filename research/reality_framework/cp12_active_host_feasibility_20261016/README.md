# GRUT: bounded full-active-host feasibility

Research base: `609db3f19576f1fe6d7713748a7564c91c3a49b7`.
Recovery checkpoint: `c2751a047c3ad57eb415b39228297b3f796b9287`.

This packet tests Path A under **unchanged frozen CP-12**, then develops
one explicitly new Path B transport proposal far enough to expose its
activation obstruction. No Path B term enters a Path A numerical run.

The full active host is a more defensible conditional transport metric
than the failed circulating ring. The five-site host-only source experiment
resolves a strong hierarchy at distances one, two and three. Original
neutral dynamics modestly bias host distributions, but a spatial phase,
continuum, gravity and a distinctive empirical prediction remain unearned.
The report includes the six-site response with its actual error intervals;
unresolved signals are never called zero or a successful asymptotic front.

## Review

- [REPORT.md](REPORT.md): generated from the completed evidence.
- [HANDOFF.md](HANDOFF.md): continuation and architecture decision.
- [PROTOCOL.md](PROTOCOL.md): bounded budget and state-independent host choice.
- [MATHEMATICS.md](MATHEMATICS.md): forced diameters, complete record work,
  exact cross-host norm identity and size-dependent persistence bound.
- [PATH_B_DESIGN.md](PATH_B_DESIGN.md): a new local flux-gated neutral-probe
  law, its formation identity and exact unsuppressed activation control.
- [REVIEW_CONTINUITY.md](REVIEW_CONTINUITY.md): checked frozen SCOUT-1 history;
  the proposed quotient-value campaign already ran and is not restarted.
- [ACTIVE_HOST.svg](ACTIVE_HOST.svg): actual graph and response data.
- [ASSESSMENT.json](ASSESSMENT.json), [VERIFICATION.json](VERIFICATION.json)
  and [manifest.json](manifest.json): claim decisions, checks and hashes.

## Reproduce

Tested environment: Python 3.12, NumPy, SciPy, Matplotlib and a C++ compiler
for the inherited signed-symmetry engine. Exact versions appear in
`VERIFICATION.json`. Run large jobs serially in about 10 GB RAM. Cold-cache
six-site matrix reconstruction takes substantially longer than checking
published outputs. Numerical caches are disposable and not committed.

From the repository root:

```sh
python research/reality_framework/cp12_active_host_feasibility_20261016/reproduce.py
python research/reality_framework/cp12_active_host_feasibility_20261016/reproduce.py --full
```

The default verifies saved evidence and the content manifest. `--full`
regenerates graph catalogs, original-neutral host marginals, host-only
source experiments, independent labelled controls, exact Path B algebra,
report and figure, then refreshes the manifest. `GRUT_SHORTCUT_CACHE`
stores inherited basis caches and `GRUT_SOURCE_CACHE` stores this packet's
neutral/source states; temporary directories are used by default.

Every unretained physical hopping and rewiring row remains in the error
action. Initial numerical-tail errors and conditional normalization are
included. These are analytic residual/integration bounds evaluated in
ordinary floating point with explicit allowances, not outward-rounded
formal certificates. Successful numerical checks do not establish any
GRUT physical milestone. Canonical theory, Stage3, PR12 and all inherited
Hamiltonian sources remain unchanged.

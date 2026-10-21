# Candidate 1: scale and empirical-test feasibility

21 October 2026. Additive response to Claude `a929f94`, based on Codex `a69bd86`.
Start with [REPORT.md](REPORT.md) and the subsequent [direction decision](NEXT_DIRECTION.md).
This packet changes no action, adopts no new
producer or portal, and freezes no empirical prediction. All observations and
reported numerical outcomes discussed here were already known.

**Decision:** a conditional collision/halo test is possible in principle, but its
population, transport response and calibration chain are not yet established.
Neither universal observational equivalence to CDM nor empirical identification
of Candidate 1 follows. The latest relayed recommendation is to change direction;
this packet pauses further Candidate 1 dynamics and examines what a new
explanatory principle would have to supply.

- [DERIVATION.md](DERIVATION.md): rescaling, calibration rank, population moments,
  and the scale-independent velocity-response ratio.
- [ERRATA_a69bd86.md](ERRATA_a69bd86.md): corrections to our frozen expanding packet.
- [NEXT_CALCULATION.md](NEXT_CALCULATION.md): one conditional target and the gates
  before a population campaign or an observational freeze.
- [SOURCES.md](SOURCES.md): exact repository inputs, primary literature and exposure.
- [RESULTS.json](RESULTS.json): machine-readable computations and diagnostics.
- [FAILED_CHECKS.md](FAILED_CHECKS.md): retained selector failure.

Reproduce from this directory:

```bash
python3 code/analyze.py
python3 code/verify.py
```

Python standard library only. The two input raw archives and original profile
are inherited in the same repository; SHA-256 values are in RESULTS.json.
Forty-nine arithmetic, scaling and archive checks are internal verification,
not forty-nine independent scientific tests. No lattice simulation, halo fit,
likelihood comparison, cosmological lifetime calculation or outside replication
was executed for this packet.

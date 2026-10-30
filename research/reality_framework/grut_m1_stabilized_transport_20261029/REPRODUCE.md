# Reproduction

Use a checkout containing this packet and the pinned parent tree. Python, NumPy,
SciPy, Matplotlib and a C11 GCC are required. Versions used are recorded in
[ENVIRONMENT.json](ENVIRONMENT.json). Run from this packet's directory:

```bash
bash reproduce.sh checks
```

This reruns the independent background, mathematical/conservation checks,
representative density response, note counterexamples and summary tables.
It reads the archived convergence runs and optional CLASS spectra. The output
checks are same-author tests, not external replication.

```bash
bash reproduce.sh full
```

This additionally rebuilds all holder/phase, elastic and internal hierarchy
cases. Single-thread BLAS avoids large overhead on the small repeated matrix
operations. Runtime depends on the host; the 512-node elastic operator and
high-multipole hierarchy are the more expensive cases. `.work/` contains
disposable compiled binaries, generated DOP853 coefficients, solver inputs
and logs; it is not required from this publication.

Optional CLASS reruns use [this exact source commit](https://github.com/lesgourg/class_public/tree/814ce59c380a5cccd0bed1747a384ff7bfc1f9cf).
Build it separately with `make class -j2`, then:

```bash
GRUT_CLASS_BINARY=/absolute/path/to/class bash reproduce.sh class
```

Both successful moderate-k controls and the high-k failure/timeout attempts
are scripted. A failure status in a high-k CLASS attempt is an archived
negative result, not an instruction to suppress the error. The successful
controls are collisionless and adiabatically back-extrapolated; they do not
validate the formation matching or BGK treatment.

Run `python code/package.py` to check pinned parent bytes, parse source and
data files, check local document links, and write the SHA-256 manifest. The
manifest excludes itself and transient files. Reproduction can change wall
times, floating-point details and archive metadata; compare numerical fields
using the declared controls rather than demanding identical regenerated
output bytes. The archived input pin hashes are exact.

File dictionary: `holder`, `reference`, `refined` and `finite_k` contain the
background/transfer data; `phase_*` are direct phase-perturbed backgrounds;
`elastic_*` carry normalized comoving energy populations, rates and fits;
`linear_*` carry conditional density histories and ratios; `warmrk*` and
`coldrk*` are successful CLASS spectra. Definitions and omitted physics are
in [DERIVATIONS.md](DERIVATIONS.md). The figure connects sampled k values;
it is not a dense scan certifying every acoustic feature.

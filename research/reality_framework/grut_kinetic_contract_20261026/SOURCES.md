# Sources and credit

The exact source records and file hashes are in [SOURCES.json](SOURCES.json).
This packet uses our published 225028b finite band; it does not add spectrum
tails or perform a new mode integration. The independent finished review is
[88fe3b0](https://github.com/ryangrvr/TestingGRUT/commit/88fe3b08d890baa0acebcb4f3fb92f487fbe3b10),
summarized on [issue #13](https://github.com/ryangrvr/TestingGRUT/issues/13#issuecomment-6409470774).
Its conditional scalar result is accepted. Source production and the full
kinetic calculation have distinct ownership; no external human replication
is implied by separately written AI integrators.

The f7d641a lsrelic.py rate and state functions were read directly from an
immutable archive. The existing scalar route has ss conversion and a
fixed-temperature inventory denominator. Its state functions omit explicit
neutral-scalar thermodynamics. The present packet imports none of that code.

Two primary papers were checked for methodological context:

- Cervantes and Hryczuk, [Freezing-in Cannibal Dark Sectors](https://arxiv.org/html/2407.12104v2),
  especially section 2 and section 5.2. They evolve coupled species moments
  and state the assumptions behind distribution closures. Their models,
  rates and results are not a validation of this GRUT action or our source.
- Ala-Mattinen, Heikinheimo, Kainulainen and Tuominen,
  [Anatomy of real intermediate state-subtraction scheme](https://arxiv.org/html/2309.16615v2),
  introduction and minimal double-counting example. It explains the overlap
  when an unstable intermediate species is tracked both as a real particle
  and through resonant scattering. We recommend a stated consistent treatment;
  no double count of the old threshold amplitudes is demonstrated here.

The decay-source transform, reaction inventories, reservoir quadrature and
finite-Q differential identity are derived or checked here. None of the
papers supplies a formation kernel, physical phase holder or empirical signal
for this branch.

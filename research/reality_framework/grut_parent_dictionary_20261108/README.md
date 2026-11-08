# P1-01 — smaller parent-cutoff changes, physical matching still unresolved

The fixed-physics extension from 1280 to 2894 states reduces every tested
successive-cutoff history difference. It does **not** meet the registered
cutoff, overflow or lattice-matching gates. The largest lattice-to-parent
discrepancy changes only from 0.00247220 to 0.00242882, against a target of 0.001.
This is a regulated comparison, not a measurement of the unknown absolute
continuum error. No physical parameter was refitted.

This packet follows the physical-bridge priority rather than adding a memory
algorithm. [DICTIONARY.md](DICTIONARY.md) derives the boundary-induced mean,
records the mass and detector conventions, and distinguishes parent truncation
from E1's additional reduction. It reviews the adopted theory; it is not an
independent external review of our code.

## Registered outcome

The equations, protocol, new wrapper and analysis were frozen publicly at
[`18ec80f0c0a608ac4edad675cf5144f168d84edd`](https://github.com/ryangrvr/TestingGRUT/commit/18ec80f0c0a608ac4edad675cf5144f168d84edd)
before the E14 Hamiltonian or histories were evaluated. The [protocol](PROTOCOL_P1_01.md)
and [machine-readable card](protocol_p1_01.json) preserve that sequence.
The old M1 results and N20 comparator were already examined; this is not an
untouched holdout. Production used a separate clean checkout of that commit.

| Quantity | Earlier M1 | P1 | Frozen criterion | P1 verdict |
|---|---:|---:|---:|---|
| Parent dimension | 1280 | 2894 | <=3000 | within cap |
| Maximum successive-cutoff TV | 0.00140056 (E10→E12) | 0.000338001 (E12→E14) | <=0.0001 | miss |
| Maximum parent overflow | 0.00126337 | 0.000969100 | <=0.00001 | miss |
| Maximum N20-to-parent TV | 0.00247220 | 0.00242882 | <=0.001 | miss |
| Quadrature refinement TV | 2.87e-14 at E10 | 2.77e-14 at E14 | <=1e-6 | pass |

The maximum successive step is 0.24133 times the preceding one, or 4.14 times
smaller, with 2.26 times the dimension. That is a convergence trend over the
tested steps, not a proven asymptotic rate or an extrapolated continuum value.
The matching discrepancy falls by only about 1.75%. It cannot yet be assigned
uniquely to the lattice, physical dictionary or unresolved bosonic truncation.

| Times | E10→E12 TV | E12→E14 TV | N20↔E14 TV | E14 overflow |
|---|---:|---:|---:|---:|
| (0) | 0.000378030 | 0.000147489 | 0.000768644 | 0.000299057 |
| (0,1) | 0.000935486 | 0.000228784 | 0.001914631 | 0.000638547 |
| (0,2) | 0.000681591 | 0.000292937 | 0.001462267 | 0.000598229 |
| (0,1,2) | 0.001400556 | 0.000338001 | 0.002428822 | 0.000969100 |

Every ordinary outcome and the absorbing overflow outcome remain in each TV
comparison. No surviving probability is renormalized. For any complete physical
distribution with no artificial overflow, TV is at least the truncated overflow
probability. At E14 that necessary floor is just below 0.001 for the longest
history, but the observed discrepancy remains above it. Overflow is not an
upper error bound and cannot be subtracted to certify physical accuracy.

## What the dictionary review resolves

1. **Nonzero boundary mean.** The declared standard Gauss law in the neutral
   even chain corresponds to phi(0)=phi(L)=sqrt(pi)/4 in the adopted continuum
   map. Solving the massless quadratic boundary-value problem yields the cosh
   profile used by the builder. It is an oscillator displacement, not an imposed
   massive expectation. The computed massive mean of F changes from 0.01789928
   at E12 to 0.01843405 at E14, while the free displacement mean is 0.12722990.
   This diagnostic has no newly invented accuracy gate and is not a fitted input.
2. **Mass and endpoint scope.** The coefficient m_R g exp(gamma)/(2 pi^(3/2))
   uses infinite-line normal ordering at mu=g/sqrt(pi). The inherited lattice
   shift m_lat=m_R-g^2 a/8 is a motivated improvement. Open boundaries break
   the translation argument's hypotheses. The L=(N+1)a endpoint placement has
   a free-hopping justification, not an exact interacting calibration. P1 does
   not silently repair either convention.
3. **Same instrument across parent cutoffs.** Padding occupation tuples gives
   an explicit nesting map. Rebuilding the E12 block inside E14 reproduces H
   within 3.55e-15 per entry and W exactly at stored precision. The positive
   decomposition of the smaller-parent overflow agrees within 2.22e-16.
   This verifies consistency of reused finite matrix elements, not an
   independent derivation of the continuum Hamiltonian.
4. **Two distinct omitted sectors.** E1 removes directions inside E12; E2
   restores those 1208 directions. P1 adds directions beyond E12. Full E2
   accuracy therefore cannot bound P1's physically separate cutoff error.

The nested-state diagnostics show 0.00230193 of the E14 ground-state norm squared
outside the E12 free-energy space; its squared overlap with the embedded E12
ground state is 0.99769360. These are changes in approximations to the same
ground-state preparation rule, not a change of the declared preparation.

For the *fixed embedded E12 ground state*, first-readout loss splits into
0.000405553 entering newly represented E12-to-E14 directions and 0.0000217512
still outside E14, summing to 0.000427304. This positive one-step identity must
not be used to allocate additive percentages of the complete history error.
The E14 ground state itself changes, as does later dynamics.

Specifically, embedding the old ground state in E14 gives first-readout overflow
2.18e-5, whereas using the consistent E14 ground state gives 2.99e-4. The high
ground-state overlap is therefore not a detector-accuracy certificate. Keeping
the old vector would understate this loss while changing the regulated
ground-state preparation rule being compared.

## Error ledger

| Contribution | Evidence from P1 | What remains uncontrolled |
|---|---|---|
| Finite matrix solver | Largest low eigenpair residual 7.43e-15; E12 replay matches all archived probabilities exactly | Independent implementation of the complete interacting model |
| Spatial integration of cosine | q128→q256 history TV <=2.77e-14 at E14 | Not a bound uniform in the energy cutoff |
| Free-mode Weyl sum | Inherited 8192-mode variance tail bound 1.91e-9; covariance quadrature estimate 4.59e-14 | These do not bound interacting preparation or dynamics |
| Parent energy/occupation truncation | Smaller successive differences; cutoff gate still misses | Absolute remainder to the interacting continuum |
| Parent instrument truncation | Explicit normalized overflow, still 96.9 times its tolerance | Omitted propagation and coherent return |
| Lattice spacing | Archived M1 N18→N20 step retained; no new lattice run | Absolute continuum extrapolation |
| Boundary/mass/detector matching | Explicit conventions and primary-literature support | Finite-spacing and boundary corrections at requested accuracy |
| Physical apparatus | Ideal smeared instantaneous instrument specified | Local causal implementation and its error budget |

All three subprocesses completed within their registered caps. Their measured
work times were 21.25 s (E12/q128), 173.33 s (E14/q128, including nesting checks),
and 279.00 s (E14/q256). Including subprocess overhead they took about 475 s,
below the 1050 s total cap. One BLAS thread was used. This is parent-level
construction and exact dense propagation cost, not an online compressed solver
or a claimed speedup.

## Reproduction and actual validation

From the repository root, install the pinned dependencies in a fresh environment:

```bash
python3 -m venv /tmp/grut-p1-venv
/tmp/grut-p1-venv/bin/python -m pip install -r research/reality_framework/grut_detector_matching_20261107/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /tmp/grut-p1-venv/bin/python research/reality_framework/grut_parent_dictionary_20261108/verify_results.py --negative-controls
```

To repeat production from the frozen source, choose paths that do not already
exist. The raw-output directory must be outside the checkout:

```bash
git worktree add --detach /tmp/grut-p1-frozen 18ec80f0c0a608ac4edad675cf5144f168d84edd
cd /tmp/grut-p1-frozen
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /tmp/grut-p1-venv/bin/python research/reality_framework/grut_parent_dictionary_20261108/run.py --output-dir /tmp/grut-p1-replay
```

`verify_results.py` is a postproduction preservation check in the result commit;
it is not present in the earlier freeze. `run.py` and `analyze.py` are frozen.
The recorded production used Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0.
Same-environment replay and stability checks are not independent replication.

Actually checked: the prefreeze small nested-block identity; production from a
clean public source; exact E12 replay; the E14 nesting/overflow identities;
all 12 complete saved distributions (largest normalization discrepancy
3.78e-15); basis/checkpoint/spectrum consistency; 13 frozen-source/input hashes;
and recalculation of every saved verdict. The verifier rejects a forged
continuum-success verdict and a history with overflow removed, even when derived
comparison arithmetic is recomputed. See [VALIDATION.json](VALIDATION.json),
[raw analysis](raw/analysis.json), [events](raw/events.json) and
[provenance](raw/provenance.json).

The bounded calculation stops here. There is no E16 extension, adaptive mass or
boundary tuning, new memory algorithm or new lattice grid in P1. E1 and E2 keep
their original finite-parent verdicts. M0/M1 misses remain on record. Physical
detector-history matching, a new GRUT restriction and observational identification
remain unestablished. No merge, release or Zenodo upload is authorized.

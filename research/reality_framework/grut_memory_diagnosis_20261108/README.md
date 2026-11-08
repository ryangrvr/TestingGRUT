# E2-01: the remaining error is mainly detector-mediated; full memory controls it

8 November 2026 · Finite-parent diagnosis and continuation · Zenodo paused

**The tested autonomous error is dominated by discarding detector-created
amplitudes outside the retained image, under the declared factorial controls.**
Restoring propagation exchange alone does not improve the old histories.
Restoring detector transfer and retaining the resulting amplitudes under their
complementary dynamics does. The signed propagation, instrument and interaction
contributions partly cancel; they are not additive percentages.

Keeping all omitted amplitudes as memory then reproduces the parent on both old
and newly frozen timings. At the registered decision resolution, its analytic
error bound is below 0.001. This succeeds by carrying **72 retained plus 1208
memory amplitudes**. It is not a compressed 72-state solution or an interacting
continuum certificate. E1-01's strong accuracy miss and both bound misses remain.

[Derivation](DERIVATION.md) · [Frozen protocol](PROTOCOL_E2_01.md) ·
[Implementation](experiment.py) · [Raw summary](raw/summary.json) ·
[Validation record](VALIDATION.json)

## Frozen scope and chronology

Protocol, equations and code were committed at
`3dce73fb6d9e6822281c048f8cc646ebdd8b74ba` before the new controls, bounds,
continuation or timings were evaluated. Only the declared small synthetic matrix
checks ran before that freeze. The production run used a clean detached checkout,
one BLAS thread, and the existing NumPy 2.3.5 / SciPy 1.17.0 environment.

The parent is unchanged: E=12, dimension 1280, m/g=0.5, Lg=5, theta=0,
Dirichlet boundary phi/sqrt(pi)=1/4 and 128 spatial quadrature nodes. E1's stored
72-column isometry is reused without a new anchor. Detector strength remains
c=2, width 0.8, phase zero. The old timings are calibration/replay; the two new
timings were frozen in this registration and were not selected from a scan.

## 1. What the controls isolate

The two switches restore dynamical exchange between the retained image and its
complement, and detector transfer between them, respectively. All four processes
are normalized instruments acting on the same prepared state. Both-off stays in
the 72-state image and reproduces E1's joint autonomous approximation. Both-on
is the original parent to the separately recorded preparation discrepancy.

| Old timing sequence | Both off: autonomous | Propagation only | Detector transfer only | Both on vs original parent |
|---|---:|---:|---:|---:|
| (0,0.7,1.9) | 0.001508985 | 0.001551782 | 0.000167176 | 3.95e-15 |
| (0,0.6,1.5,2.4) | 0.002290737 | 0.002393296 | 0.000230135 | 6.11e-15 |

Entries are total-variation distances, not relative errors. “Detector transfer
only” retains the entire complementary state and its C evolution; it does not
mean an autonomous 72-state detector correction. It includes later detector
return from that state. The counterfactuals change declared operations and are
diagnostic models, not four equivalent physical descriptions.

For the longer old history, the half-l1 norms of the signed changes from the
both-off model are:

| Signed component | Half-l1 norm |
|---|---:|
| Propagation exchange | 0.000398555 |
| Detector transfer | 0.002312936 |
| Interaction of the two switches | 0.000279346 |
| Norm of their signed sum | 0.002290737 |

The component norms sum to 0.002990837, larger than the norm of their sum.
Every signed outcome difference, including overflow, is archived and the
factorial identity reconstructs exactly at stored precision. Thus “detector
transfer is dominant in these controls” is supported; a universal percentage of
the error caused by each mechanism is not.

Preparation is negligible here: the embedded state's discrepancy from the parent
ground vector is 2.46e-15 and its eigen-residual is 7.30e-15. This does not remove
the detector's ability to populate excited complementary states.

## 2. Which branches matter

The first detector sends a total ordinary-branch norm squared of
0.000519689 into the complementary image, despite the negligible initial
dynamical residual. Subsequent instruments and evolution can return those
amplitudes. An autonomous truncation instead adds their loss to its overflow.

The largest unnormalized autonomous detector-leak norm is 0.0169275 on prefix
`++`, immediately before the t=1.5 readout of the longer old history. This is an
amplitude norm, not a probability. That prefix also has the largest dynamical
source norm, 0.061886 in the declared units, and weight 0.932184. Lower-weight
prefixes can have larger conditional residuals; for example, `+-` at t=1.9 has
weight 0.03656 and conditional residual 0.24113. No such branches are discarded.
[All rankings](raw/rankings.json) and parent prefix diagnostics are retained.

For the longer old history, autonomous overflow is 0.00411582 versus the parent's
0.00183825. Much of the total history discrepancy appears as this extra loss.
The parent's own overflow remains; the continuation does not solve the parent's
missing continuum modes.

An explicitly normalized control erases complementary amplitudes after every
nonfinal readout and transfers their norm into absorbing overflow. Its maximum
TV is 0.00228137 on old timings and 0.00530607 on new timings. Restoring exact
propagation while repeatedly deleting the detector-created memory therefore
does not meet the 0.001 target in either group. This reset is a changed instrument,
not an innocuous numerical operation.

## 3. The continuation and its finite-parent bound

The memory equation in the derivation is implemented through its equivalent
coupled p,q amplitudes. It preserves q across all readouts. Block-diagonal
evolution and retained/complement exchange are applied by unitary Strang steps;
the exchange uses every singular value of the coupling matrix. No projected
detector is made artificially unitary. The full parent detector keeps its positive
overflow effect and complete instrument.

| Timing group | Autonomous maximum TV | Memory maximum TV at n=512 | Maximum analytic bound at n=512 |
|---|---:|---:|---:|
| Old calibration | 0.002290737 | 3.46e-9 | 0.000067998 |
| Newly frozen timings | 0.004174062 | 5.71e-9 | 0.000092482 |

The new timings are (0,0.43,1.37,2.61) and (0,0.31,1.07,1.83,2.79). Both groups
pass the registered measured-accuracy, method-improvement and analytic-bound
thresholds for this memory continuation. They are new numerical configurations,
not independent observational data. The bound is conditional on the finite
supplied model and is evaluated in ordinary floating point. It is not an
interval-arithmetic certificate.

All fixed refinements are reported:

| Substeps per nonzero interval | Old maximum TV | New maximum TV | Largest bound across all histories |
|---|---:|---:|---:|
| 32 | 8.87e-7 | 1.46e-6 | 0.0236755 |
| 128 | 5.53e-8 | 9.14e-8 | 0.00147972 |
| 512 | 3.46e-9 | 5.71e-9 | 0.0000924822 |

The measured errors decrease approximately sixteenfold when the number of steps
is quadrupled, consistent with second-order global error. Lower-resolution bound
misses are retained. The decision level n=512 was fixed before execution.

The bound uses the established two-term real-time Strang estimate from
[Childs et al., equations (146)–(150)](https://yuansu.me/publication/trotter_error/trotter_error.pdf).
The two Frobenius commutator upper bounds are 157.26099 and 24.65434, giving
c_S=8.6070694. Complete-instrument dilations and trace-distance contraction
extend the unitary error bound to the full history distribution. This is standard
quantum simulation mathematics applied to the adopted model, not a new GRUT law.

## 4. The tighter autonomous bound still fails

The new branch-relevant bound uses exact parent propagation defects on every
autonomous branch and combines branch norms before bounding probabilities.
It keeps preparation, propagation and instrument terms explicit and omits none.

| Timings | Preparation norm | Accumulated propagation defect | Accumulated instrument defect | TV upper bound |
|---|---:|---:|---:|---:|
| (0,0.7,1.9) | 2.46e-15 | 0.036277 | 0.070017 | 0.212435 |
| (0,0.6,1.5,2.4) | 2.46e-15 | 0.066207 | 0.099712 | 0.331497 |
| (0,0.43,1.37,2.61) | 2.46e-15 | 0.070524 | 0.110080 | 0.360782 |
| (0,0.31,1.07,1.83,2.79) | 2.46e-15 | 0.115108 | 0.148416 | 0.526187 |

These positive budgets are upper-bound ingredients, not the signed attribution
in section 1. They dominate the observed errors but do not establish 0.001
autonomous accuracy. Their exact-propagator evaluation also incurs parent cost.
The successful bound in section 3 belongs to the different memory-retaining
approximation and cannot be transferred to E1.

## 5. Verification, cost and scope

The bounded production calculation took 29.4 seconds on one BLAS thread. It
reused the E1 isometry and parent builder, constructed a complementary basis,
diagonalized its C block and the parent reference, and evolved every auxiliary
amplitude. Full-parent eigenvectors are used for exact controls and diagnostic
bounds, not to propagate the low-rank split step. Dense parent matrices and
full auxiliary state storage are still required. No computational speedup or
minimal sufficient memory has been established.

Checks actually performed:

- Complex small-matrix exchange versus independent direct exponentiation,
  commuting limit, complete overflow dilation and local commutator bound.
  An expansive detector is rejected.
- Parent/autonomous E1 replay, coordinate isometry, Hamiltonian reconstruction,
  Hermiticity, instrument contraction/completion and the compressed-overflow
  identity. Coordinate residual is 1.92e-13; Hamiltonian reconstruction residual
  is 1.09e-12. These are numerical diagnostics, not interval enclosures.
- Recomputed all 32 stored distributions' normalization and TV, all signed
  factorial terms, bound arithmetic, source hashes and frozen gate decisions.
  Largest distribution normalization error is 4.89e-14. A deliberately corrupted
  zero bound is rejected by the separate output verifier.

These checks do not constitute an independent implementation of the full physical
calculation. Output integrity and reproducibility also do not establish that the
research evaluator reliably judges unfamiliar theories. Those remain separate
questions. The boundary/mass-dictionary review has not been completed by this test.

## Reproduction

Use a clean repository checkout containing this packet, Python 3.12 with
NumPy 2.3.5 and SciPy 1.17.0. From the repository root:

```sh
export GRUT_MEMORY_ROOT="$PWD"
export GRUT_MEMORY_OUTPUT="$(mktemp -d)"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
timeout 600 python "$GRUT_MEMORY_ROOT/research/reality_framework/grut_memory_diagnosis_20261108/experiment.py" --output-dir "$GRUT_MEMORY_OUTPUT/e2"
python "$GRUT_MEMORY_ROOT/research/reality_framework/grut_memory_diagnosis_20261108/verify_results.py"
```

The first command refuses a dirty checkout or existing output directory. For
exact freeze replay, check out `3dce73fb6d9e6822281c048f8cc646ebdd8b74ba` and use
the first command; the later report/verifier need not be present. The second
command checks the committed raw output by default; it can inspect regenerated
output with `--raw "$GRUT_MEMORY_OUTPUT/e2"`. It requires the scientific source
hashes to equal the frozen run's hashes and reports whether the checkout commit
itself is the freeze or a later documentation/results commit.

A completed run may contain failed scientific gates: read `summary.json`.
Exceptions are preserved as `FAILURE.json`; timeouts preserve stage outputs and
must not be called passes. Nothing is merged or released. **Zenodo remains paused.**

## Result retained for GRUT

This finite realization now identifies a concrete failure of an autonomous
description: detector-generated complementary amplitudes retain predictive
information across readouts. The explicitly retained memory restores the tested
histories with a finite-model error bound, at full state-space cost. The
construction makes preparation, dynamics and instruments consistent without
deriving their physical selection. Continuum matching and distinctive GRUT
physics remain open, with the earlier failures preserved.

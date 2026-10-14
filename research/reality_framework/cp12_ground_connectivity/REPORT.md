# CP-12 — actual ground connectivity, graph mixing, and an electric-tail certificate

14 October 2026 UTC. Work-led additive research on TestingGRUT, based on
connected-bridge checkpoint 3b914e5a730a8d246175159cf5bd4d3f6870404a.
Frozen parameters, Q=N, zero-temperature ground states.

## Finding

**The small full quantum ground state is predominantly connected. A
disconnected variational trial does not establish a fragmented ground state.**
The calculation also identifies why this is insufficient for space:
connected graphs with short cycles and attached trees carry much more
probability than a single extended ring.

This pass includes negative fermion hopping, Gauss-valid integer fields,
the rewiring diagonal, and coherent transitions between all graph shapes.
The six-event result is checked in all eleven permutation representations,
covering a regulated physical Hilbert space of **36,801,265 states**.
It is not obtained by assuming the ground state is a permutation singlet.

## Ground-state results in one declared regime

All entries below use Q=N and electric cutoff |E|<=3.

| N | Lowest energy | Connected probability | Uniform-graph baseline | Single spanning ring | Mean fraction in largest component |
|---:|---:|---:|---:|---:|---:|
| 4 | -0.598552282476 | 100% | 100% | 28.9561% | 100% |
| 5 | -0.728493363556 | 94.9037% | 88.0952% | 9.48709% | 98.9807% |
| 6 | -0.858318930157 | 86.6971% | 73.1269% | 3.22019% | 97.2769% |

The baseline uniformly counts labelled graphs with N edges; it is not a
CP-12 thermal ensemble and does not weight graphs by their matter/rotor
Hilbert-space dimension. At N=4 every permitted graph is connected, so
its 100% is a combinatorial fact, not a selection result.

The ground has enhanced connectedness relative to this baseline. Both the
enhancement and declining absolute probability are finite-system facts.
Neither establishes a thermodynamic limit.

At N=6, the probability that every degree equals two is 3.76708%. This
includes disconnected pairs of triangles; the probability of one ring is
only 3.22019%. Conditional on connectedness, the sole cycle length has
distribution:

| Cycle length | Conditional ground probability |
|---:|---:|
| 3 | 50.7848% |
| 4 | 31.2796% |
| 5 | 14.2213% |
| 6 | 3.71430% |

The mean cycle length conditional on connectedness is 3.70865. The finite
ground predominantly joins vertices through short cycles with attached
trees. We have not derived an extended spatial phase.

## Checking the missing symmetry sectors

The previous invariant solver gives a candidate ground in one symmetry
sector. The new solver retains all signed state orbits, constructs every
Young representation of S_N, projects its stabilizers, removes null
projector directions, and solves every multiplicity block.

Independent controls:

- Complete N=3,L=2 spectrum: 78 eigenvalues reproduced within 1.3e-14.
- Complete N=4,L=1 spectrum: 1,875 eigenvalues reproduced within 8.5e-14.
- Direct N=5,L=3 full-matrix solve: 427,986 states; agrees with the
  representation reduction on ground energy, next energy, and connectedness.
- N=6,L=3: all eleven sectors solved; dimension-weighted multiplicities
  sum to exactly 36,801,265.

The six-event invariant ground energy is -0.858318930157. The lowest
nontrivial representation, [4,2], has energy -0.787204005610. The invariant
sector's next energy is -0.686069023458. The resulting numerical full gap
is 0.071114924547. This gap identifies the finite ground; it is not a
momentum dispersion or evidence of a graviton.

Large eigensolves use floating arithmetic and checked residuals. They are
numerical spectral identifications, not exact rational proofs that every
omitted Lanczos direction lies higher.

## A projector bound that includes graph mixing

Let E_C and E_D be connected/disconnected block minima, b the norm of their
cross coupling, and E the ground energy in the declared symmetry sector.
For delta_C=E_C-E>0 and delta_D=E_D-E>0, the exact eigen-equations imply

\[
\boxed{\frac{\delta_D^2}{\delta_D^2+b^2}
\le P(C)\le
\frac{b^2}{\delta_C^2+b^2}.}
\]

For the actual invariant ground components, numerical evaluation gives:

| N | Lower comparison | Computed P(C) | Upper comparison | Cross-sector mixing energy |
|---:|---:|---:|---:|---:|
| 5 | 90.2560% | 94.9037% | 97.5473% | -0.0248697943 |
| 6 | 77.4165% | 86.6971% | 92.5753% | -0.0610164317 |

Both block matrices retain the original rewiring boundary diagonal. The
six-event block minima quoted in the data are invariant-sector minima;
projecting its ground onto connectedness preserves that sector. Their
use in the bound does not assume that every restricted block minimum is
also its full-space minimum.

The negative mixing energy is a direct contribution to the ground energy.
It is lost if one dephases only the binary connected/disconnected
information. For N=6 that mathematical operation raises the system
energy by 0.0610164317. No physical measurement apparatus is supplied.

This resolves the immediate logical gap in the bridge audit: its positive
bridge cost alone cannot determine the complete ground-state competition.
The negative hopping and coherent coupling actually matter.

## Electric cutoff: numerical agreement and an exact small certificate

At N=4,5,6, increasing L from 3 to 4 in the invariant block changes the
ground energy and connectedness by less than 1e-10. Boundary probabilities
fall from about 1e-11 at L=3 to about 1e-19 at L=4. This is a numerical
regulator check; it is not an infinite-cutoff proof. All representations
are checked at L=3, not at L=4.

We also derive a route to a proof. Let A_L be the retained block, C_L the
omitted block, and B_L the coupling from retained to omitted states.
The frozen star-energy bound and paired-leg Gauss law imply

\[
C_L\succeq
\left[N(\tfrac32-\sqrt3)+\frac{(L+1)^2+L^2}{8}\right]I=c_LI.
\]

An outside state reachable in one hop has a unique inside predecessor.
Therefore \(B_L^\dagger B_L=D_L\) is a diagonal outward-hop counter.
For z<c_L,

\[
A_L-zI-\frac{D_L}{c_L-z}\succeq0\quad\Longrightarrow\quad H\succeq zI.
\]

At N=3 this comparison is verified by exact rational elimination. Together
with a rational trial-state upper bound it proves, for unbounded integer
electric fields,

\[
\boxed{-0.46729956\le E_0\le
-\frac{18691982041980411}{39999999816851996}
\approx-0.4672995531891348.}
\]

The interval width is 6.81e-9. This is an actual small unregulated
certificate, not merely agreement between two electric cutoffs.

The same comparator was tested numerically in every representation:

| N | Threshold z | Smallest comparator eigenvalue | Sectors checked |
|---:|---:|---:|---:|
| 5 | -0.728494 | +6.36444e-7 | 7 of 7 |
| 6 | -0.85832 | +1.06984e-6 | 11 of 11 |

These positive numerical minima support the cutoff interpretation, but
the large matrices have not received exact or interval-certified
positivity tests. Their result must not be promoted to the status of the
N=3 rational certificate.

The Schur technique is established mathematics. The CP-12 application and
diagonal outward-hop identity are derived in DERIVATION.md. See Dusson,
Sigal, and Stamm, [arXiv:2105.02058](https://arxiv.org/abs/2105.02058) for the
general Feshbach-Schur framework.

## Provisional larger checkpoint

An N=7,L=2 invariant solve is retained in the raw data: connected
probability 77.4268%, single ring 1.11379%, and mean largest-component
fraction 95.5840%. Its boundary probability is 2.50e-5. Other
permutation sectors and cutoff refinement were not computed. It is
excluded from the validated full-ground table and from large-system
extrapolation.

## Scientific decision

The claim that the actual CP-12 ground must fragment is not established
by the disconnected energy witness, degree-two exclusion, or bridge
penalty. The finite full-ground calculation instead finds a mostly
connected but predominantly non-ring state.

The same data warn against treating exact connectedness as the whole
question: a state can fail to be completely connected while retaining
most vertices in one component. The next large-system test must track
that component's intrinsic volume growth, cycle structure, and physical
propagation together in the same state. A fit to four connectedness
percentages would not settle it.

This packet adds a complete symmetry-sector calculation, a graph-projector
comparison, and a small unbounded-rotor certificate. It does not establish
a thermodynamic geometry, a spacetime continuum, or Einstein gravity.

## Reproducibility

See README.md for commands and VERIFIED.json for exactly which checks ran.
The default verifier recomputes small spectra, exact certificates, and
graph counts, then checks saved large results and their cross-agreement.
It does not silently rerun every large eigensolve. Large solvers have
an explicit --fresh option to bypass result resumption.

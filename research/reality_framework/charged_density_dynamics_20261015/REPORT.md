# GRUT builder: charged-network formation, return and persistence

This executes the analytical continuation at `1e68b1a`, after the memory
bridge `5d63406`. The frozen CP-12 law and original neutral preparation are
unchanged. The main new result is finite-time charged circulation density
with controlled outgoing-state errors at N5 and N6, supplemented by exact
local onset kernels through N8. Spatial phase selection is still open.

## The definition repair matters

The proposed equality E_ie=-E_je excluded every edge hole and made its
all-edge charged event empty. The repaired primary criterion is opposite
nonzero signs, as in the published charged-circulation prototype. A hole may
carry (2,-1). Equal-and-opposite circulation is tracked separately and obeys
c_bal<=Q-n_v. The amendment to the incoming note is explicit.

The declared primary event requires circulation on at least one quarter of
all active edges and vertex matter on at least one quarter of vertices.
No degree, host topology or spatial dimension is selected. At N5,N6 this is
c>=3,n_v>=2. It is a weaker event than the old all-edge four-regular target;
the new probabilities do not replace that target's much rarer probabilities.

## Finite-time formation now has positive lower estimates

| System and time | Central probability | Full-law probability interval |
|---|---:|---:|
| N5,Q10,T1 | 0.00264308339 | [0.00264308097, 0.00264308582] |
| N5,Q10,T2 | 0.05962997288 | [0.05961615929, 0.05964378808] |
| N5,Q10,T3 | 0.26794386049 | [0.26577409877, 0.27012244314] |
| N5,Q10,T5 | 0.46103962729 | [0.39893044726, 0.52764035167] |
| N6,Q12,T0.5 | 0.00000718418 | [0.00000608258, 0.00000837742] |
| N6,Q12,T1 | 0.00357020302 | [0.00040949290, 0.00985382229] |

The N5 electric-square budget is 50, with 162,191 invariant coordinates.
The N6 budget is 12, with 715,985 invariant coordinates representing
513,177,325 physical orbit states. Every outgoing unregulated transition
enters the error estimate. This is not a 513-million-state diagonalization.
The estimates evaluate analytical Duhamel bounds in floating point with a
separate reporting allowance; they are not formal interval-rounding certificates.

N5 reaches appreciable partial circulation. Its stricter half-edge threshold
also has probability 0.25644 at T5, with interval [0.21069,0.30668]. Requiring
three quarters gives [0.00672,0.03125]. The full-edge T5 readout is unresolved.
All sensitivity thresholds and the broad T10 estimate are retained. The latter
has interval [0.04135,1], not a precise long-time probability or current.

N6 has genuine graph choice: all 455 labelled Q12 graphs are in its input.
Its lower budget-eight pilot fails the formation error gate even though norm
conservation is excellent. Raising the budget resolves a positive event
weight through T1. At T2 the budget-twelve interval still includes zero; that
later-time calculation remains inconclusive.

These are fixed-fraction event definitions, but both N5 and N6 need only a
triangle and two charges. Two such small points do not establish asymptotic
extensive circulation. K5 is the only N5 host, and all N6 hosts have diameter
two. Neither is a dimension-selection experiment.

## Formation, return and history interference are computed together

At N5,T2 the event probability current decomposes as

| Contribution | Value |
|---|---:|
| Initial-source current | +0.16040039 |
| Coherent-return current | -0.07441255 |
| Total current | +0.08598784 |

The return is resolved under the action-based uncertainty bounds. This is
an exact unitary source/return decomposition, not a fitted irreversible rate.
Rewiring commutes with the circulation/matter event, but affects source and
return through the evolution inside both event blocks. Independent central
differences check the derivative convention.

At N5, wait2 and lag0.25, sequential charged-event survival is
0.73514 with conservative interval [0.73416,0.73612]. The unmeasured later
event probability is 0.08429069; recording the initial event flag gives
0.06952857. Their difference, 0.01476213, is a resolved history-interference
contribution. The normalized history overlap is about 0.322.

At N6, wait0.5 and lag0.25, the sequential joint event has interval
[1.009e-6,1.288e-5], giving a survival interval [0.1204,1]. Its raw survival
is 0.734 and raw normalized coherence 0.664, but the history-overlap error
does not resolve a nonzero interference sign there. Those central numbers
must not be promoted to controlled persistence precision or decoherence claims.

Survival concerns the charged-density event, not an identified geometry.
An initial projective record changes the later experiment. The unmeasured
probability retains both histories and their interference.

## An exact local calculation reaches larger onset targets

The small-scale onset is now derived from signed local amplitudes and exact
host counts, without constructing the N7/N8 evolution spaces:

| N,Q | Leading probability | Exact coefficient |
|---|---:|---:|
| 5,10 | T^16 | 10264/6615 |
| 6,12 | T^16 | 2228176/601965 |
| 7,14 | T^20 | 8852523251/8726130000 |
| 8,16 | T^20 | 4312328/2095875 |

The first two match independent integer orbit recurrences and complete-matrix
integer action. The local calculation reproduces their entire squared norms,
then extends to the square-loop minimum required at N7/N8. All compatible
hosts are counted; none is supplied as an initial spatial graph. The full
four-case local calculation takes under one second here. These are exact
short-time coefficients with the next even-order remainder, not probabilities
at T1 or a macroscopic phase prediction.

## Records and matter response

The complete graph/event record is evaluated without expanding every labelled
branch out of the invariant quotient. It reproduces the older N4 calculation:
energy change 0.065473972747, system entropy 4.082400908305 bits. Differences
from its full physical-basis implementation are below 8e-15 and 3e-13.

For N5 at wait2, the event record changes system energy by
0.1127110 +/-0.0011315 and produces about 0.3259743 bits of measured-system
entropy. For N6 at wait0.5, the complete graph/event record changes energy
by 0.1674930 +/-0.010481 and gives 8.829856 +/-0.009616 bits. These account
for the diagnostic intervention; they do not derive an autonomous apparatus
or a reset cost. The energy uncertainty uses the bounded change operator,
not an invalid bound on unbounded electric energy.

After an explicit graph record, the generic single-host theorem gives
same-graph survival >=0.9999 for the N6 lag0.25 interval. Random hosts have
the same normalization protection, so this is not a spatial preference.

The separate fixed-receiver response control is stored in `RESPONSE_N5_B50.json`.
It uses the same original neutral preparation and the same initial charged
event at wait2, applies exp(-i theta n0), and measures a fixed other vertex's
density after lag0.25. It retains the event probability and outgoing errors.
The unnormalized joint susceptibility is 5.72625e-5 with interval
[1.85476e-5,9.59774e-5]. Conditional on the initial charged event its value is
0.00096030, with interval [0.00031097,0.00160993]. Both exclude zero. Independent
finite source pulses agree with the commutator calculation to about 1e-11.
The receiver does not follow a final graph-distance bin. K5 still has no
distance hierarchy, so this is a matter-response control, not a causal cone.

## What this advances and what is still needed

The builder now has an executed, dimension-neutral finite-time construction
joining charged formation, coherent return, two-time event persistence,
costed records and matter diagnostics. It also has a tractable exact local
onset calculation at the next sizes. The primary finite-time formation gate
has positive lower estimates at both tested sizes; appreciable extended
autonomous spatial geometry has not been established.

The next useful numerical reduction should expand the charged boundary using
its actual outgoing residual, rather than close an entire larger electric
budget. The N6 conditional history is limited by omitted electric states,
not by sparse-exponential norm drift. Retaining destination identities and
building a charged-port reservoir is the concrete next implementation target.

The next physical test must compare a recorded/protected core's intrinsic
separation with matter response on those same states, including possible
shortcuts through the other active connections. Neutral completions cannot
be silently removed because the microscopic hopping acts on them too.

Quantum mechanics, the Born rule and external Hamiltonian time remain inputs.
No Lorentz limit, Einstein dynamics, Standard Model spectrum, universal scale
or experimentally calibrated deviation has been derived. CP-12 remains a
candidate within GRUT; successful construction, rather than preserving this
candidate, is the objective. The packet is ready for GRUT-RAI review, not
externally certified by it.

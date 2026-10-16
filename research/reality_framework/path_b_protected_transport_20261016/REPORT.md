# GRUT Path B: protected transport needs protected geometry

Base: `678e934367ead8e03cb51e39dd078ca3dc265cae`. Frozen CP-12 and
SCOUT-1 remain unchanged. GRUT-WM1 is an explicitly new proposal.

The fixed-geometry transport module now has a complete, size-uniform
virtual-path bound, a derived spectral gap, localized exact dressing
and strictly nonzero nearby transport. **The full coupled WM1 law
has not earned those properties or a spatial phase.** A separate
gapped geometry-change control proves why its gap is insufficient.

## Conditional construction: all orders, with a derived band

For a supplied static graph G of maximum degree at most four, use
`H_G = [[0,-g I],[-g I,Delta I-v A_G]]`, with Delta=8,v=1,g=1/4.
The original microscopic L basis and the exact dressed low band
are kept distinct. The exact bare Schur complement includes every
mediator walk and the dressed Hamiltonian is energy independent.

| Quantity | Analytic size-independent guarantee |
|---|---|
| Spectral gap | >= 4 |
| Exact dressed coupling | <= 0.03175416345 exp(-d/2.466303462) |
| Every nearest dressed coupling | >= 0.0009150457984 |
| Bare high-register probability, every time | <= 0.015625 |

These statements follow analytically for the declared graph family.
Numerical checks cover six cycles through N=256 and five supplied
square tori through N=256, with full microscopic source controls
at N=16,32,64. All source variances exceed 0.24 and +/-pi/2
pulse identities agree within the declared numerical allowance.
No source packet or metric in these controls is claimed to emerge.

## Stronger negative control: virtual geometry, not gap closure

GC1 is a separate four-band toy geometry, fully defined in MODEL.md.
The initial geometry is a cycle; a detuned virtual geometry has
one antipodal chord. Both graphs have degree at most three.
All original matrix elements and all intermediate states are kept.

* Full spectral gap >= 4.5, at every N.
* Starting in the bare reference band, original geometry probability >= 97.590361446% at every time.
* For every allowed low energy, antipodal Schur coupling >= 9.43805963776e-07, independently of N.

The last statement is an **all-orders lower bound**, not a leading
Taylor coefficient. A positive inverse-walk expansion makes all
additional virtual paths reinforce the four exhibited paths.
Distance is the initial cycle distance, declared before the solve.
It grows without bound while the coupling has a positive floor.
Thus this model fails the requested exponential criterion in that
metric despite its uniform gap and high geometry survival.

| N | Antipodal distance | Full Schur coupling at E=0 | Nearby coupling |
|---:|---:|---:|---:|
| 8 | 4 | 5.2492927074e-06 | 0.00102798345473 |
| 16 | 8 | 1.0245492484e-06 | 0.00102797683651 |
| 32 | 16 | 1.02344489207e-06 | 0.00102797683593 |
| 64 | 32 | 1.023444892e-06 | 0.00102797683593 |
| 128 | 64 | 1.023444892e-06 | 0.00102797683593 |
| 256 | 128 | 1.023444892e-06 | 0.00102797683593 |
| 512 | 256 | 1.023444892e-06 | 0.00102797683593 |

The exact canonical dressed-band calculation through N=64 also
retains the distant term (about 1.015e-6 at N=64); this numerical
cross-check does not extend the analytic Schur lower bound to a
separate all-size dressed-band theorem. Twelve robustness controls
vary the geometry coupling and detunings. At zero toggle coupling
the chord is inaccessible. Every positive tested toggle retains
the analytically certified nondecaying Schur contribution.

A second exact control normalizes all complete-graph hopping by
1/(N-1). Its uniform gap and row budget give only a 1/N distant
tail in an independently supplied ring metric, not exponential decay.

## Actual WM1 extension of CP-12: a surviving activation word

WM1 adds neutral low and mediator modes to CP-12 and couples the
mediators through the local flux gate. It preserves the entire
original neutral formation trajectory in probe vacuum exactly.
It also retains a five-step process that opens a charged chord,
transports the probe and restores the original geometry:
`<s,L_j|H_WM1^5|s,L_i> = -4 g^2 v`.

All 30 operator orderings were checked in exact rational arithmetic
at N=8,12,16,24 and ring flux 1,3,8, using the unregulated frozen
CP-12 transition operator. The bare-band Feshbach word itself is
nonzero. Independent exact full-H meet-in-the-middle calculations
confirm the coefficient and all four vanishing lower powers,
including every diagonal and rewiring term. These are
supplied-state, fixed-N Taylor statements.
They are **not** an all-orders large-N no-go for actual CP-12 or
a new original-neutral receiver experiment. GC1 is explicitly a
different control; its stronger theorem is not transferred to CP-12.

## Formation and acceptance decision

| Gate | Result |
|---|---|
| Explicit separate WM1 Hamiltonian | COMPLETE |
| All-orders protected static module | PROVED, supplied geometry |
| Geometry-change control satisfies desired locality | FAIL, rigorously in GC1 |
| WM1 preserves original neutral formation in probe vacuum | EXACT IDENTITY, old probabilities unchanged |
| WM1 protected coupled spatial phase | NOT ESTABLISHED |
| Formation of dressed probes from original neutral preparation | NOT ESTABLISHED; injection changes sector |
| Universal low-energy matter locality | NOT ESTABLISHED; original CP-12 a,b remain |
| Distinctive empirical prediction or earned quotient value | NONE |

A spectral gap cannot by itself increase initial low-energy
spectral weight, which is conserved by an autonomous Hamiltonian.
A viable formation mechanism must fit the preparation, energy
budget and any explicitly modeled reservoir. Merely adding a
probe detuning has not solved that problem.

The next architecture must control the **complete virtual geometry
sector**: prove that its weighted path budget or connected transport
expansion decays in an independently defined metric. Shortcuts must
require increasingly long virtual processes, or have a proved
selection/cancellation mechanism. A constant local activation
barrier is insufficient in the demonstrated control class.

## Reproducibility and limits

RESULTS.json contains all raw controls. MATHEMATICS.md supplies
the proofs, MODEL.md the input ledger and Hamiltonians, and
HANDOFF.md the next construction gate. The derived inequalities
are analytic; matrix checks are ordinary floating point with
declared tolerances, not outward-rounded certificates or external
peer review. The near-zero dressed entries of large static graphs
are numerical noise floors, never reported as detected signals.
Standard effective-Hamiltonian methods are credited in SOURCES.md.
No dimension, continuum, relativistic limit or gravity is derived.

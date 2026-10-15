# Independent in-sprint mathematical review

15 October 2026 UTC. This review was performed by a separate research worker
within the same sprint. It is not external peer review or interval-arithmetic
certification. It independently examines the new construction; it does not
certify every prior GRUT result.

## Verdict

The strongest supported construction is a well-defined quantum graph–gauge
model with three newly connected ingredients: exact charged-network formation
amplitudes from the original neutral ray, internally generated virtual gauge-loop
vertices, and a rigorous conditional bridge from a persistent graph record to
local matter dynamics under the full frozen Hamiltonian. These are substantial
constructive results. They do not yet establish that one state produced from the
neutral ray carries appreciable, persistent, three-dimensional local geometry.

No mathematical reason was found in this review to change the frozen CP-12 law.
A successor introduced solely to install a lattice or plaquette preference would
supply the result the project is trying to derive. The useful next construction
is the finite-time charged effective description, retaining actual state weights,
vertex matter and graph-motion residuals together.

## Charged-access calculation

The shortest-hop reduction is valid. Electric L1 changes by at most one per
matter hop; every other frozen term preserves it. In the minimum L1 target with
two vertex fermions, every one of d=2Q+2 events must therefore be a hop. The graph
cannot change in this term. The amplitude comes from that target graph's unique
neutral parent, with original coherent amplitude 1/sqrt(binomial(M,Q)). It is
not obtained from an unrestricted spectral norm or a different prepared state.

The Slater row-replacement identity follows by factoring the edge-to-edge block
and taking the relevant minor of its inverse product with the vertex block.
Sylvester's determinant identity supplies the determinant on the N-vertex block,
including when edge and vertex dimensions differ. Its use near T=0 is legitimate
as a formal series. At the minimum field order, diagonal degree insertions spend
hops without producing required edge variables, so only the oriented off-diagonal
variables contribute. Extracting amplitudes before taking their modulus square
is essential and is done correctly.

The cover enumeration correctly permits repeated vertices but not repeated
edges. A closed directed trail with l distinct edges has l cyclic roots even if
vertices repeat, giving l*w_l. Marker cycles supply the endpoint-permutation sign
of the principal minor. The local port pairings enumerate unordered trail/path
families with the exponential's combinatorics. This is a signed amplitude sum,
not a path-count approximation.

The factor four for two holes is justified. On a degree-four minimum-field target,
the baseline electric divergence at each vertex is even. Gauss law then forces
each of the two occupied vertex modes to receive exactly one additional positive
half-leg increment; the holes cannot share a positive endpoint. Each of those
vertices has two positive outgoing legs. Their two independent choices give four
physical targets of equal magnitude, while Fock row-order signs may differ.

The normalization and graph classifications are correct. For N6 the complement
of a four-regular graph is a perfect matching, hence 15 labelled octahedral hosts
among binomial(15,12)=455 neutral parents. For N7 the complement is two-regular,
so it is C7 or C3 disjoint C4, with 360 and 105 labelled hosts. Both complements
are connected. Their automorphism factors agree with N!/|Aut|. These two classes
exhaust the finite target, without implying that they approximate dimensions.

The positive coefficients prove nonzero target overlap with the original ray's
cyclic subspace. They do not compute appreciable occupied spectral weight,
long-time target occupation, a finite-time lower bound at growing N, or dimension
selection. The stated qualification is necessary. Leading-order inter-graph
coherence is absent because no graph changes in the shortest process; a mixture
of the same neutral parents has the same leading probability. Internal fermion
interference remains essential.

### Additional independent calculation

`algebra_check.py` numerically checks the Slater identity for three independent
complex-angle examples with one, two and three holes. Absolute magnitudes agree
to at most 3.1e-16. It also implements an independent square-free polynomial ring,
expands the cosh determinant and tanh principal minors directly, and compares
with the port-cover coefficients on N3 and N4 directed cycles with zero, one and
two marks. All eight exact rational results agree. The tanh coefficients are
independently derived by formal series division rather than imported from the
cover routine. These are independent small algebra checks, not a repetition of
the full N5/N6/N7 orbit sums or a full Hamiltonian evolution. Saved output is
`ALGEBRA_CHECK.json`.

## Conditional persistence and locality

The local hopping bound ||K_e||<=sqrt(2)|t| follows from a normalized bright
fermion mode with rotor-valued commuting phases. The norm may be evaluated in
the unconstrained space and then restricted to the invariant physical space.
Unbounded electric energy does not make this hopping unbounded.

The field-distance Dyson estimate needs no exponential prefactor. Expand around
H-K_e up to d-1 events and use the exact unitary evolution on the final remainder;
the ordered integration simplex has volume |T|^d/d!. All other terms preserve
that edge's half-leg fields, including neutral rewiring at zero field. This
justifies the final-readout protection bound and its resource cost. It is not a
continuously monitored history statement.

The fixed-host outgoing Gram identity is correct: for one initial graph g,
different pairs (removed e, inserted f) yield different destination graphs, so
the corresponding ranges are orthogonal. Each allowed replacement is an
isometry on its removable-edge subspace. Consequently

    ||(1-P_g) H P_g v||^2
      = (Gamma/M)^2 (M-Q) sum_e ||r_e v||^2.

The diagonal escape potential a*sum_e r_e, a=Gamma(M-Q)/M, must remain in H_g.
It is included. Duhamel then gives the state and response comparisons, including
C_N=Gamma*sqrt(Q(M-Q))/M and its O(N^(-1/2)) scaling at Q proportional to N.
This gain is valid for one host and cannot be averaged across an unrestricted
coherent graph superposition: outgoing amplitudes of different hosts can merge.
The manuscript explicitly preserves that restriction.

The local compression is local on the incidence graph: vertex masses are onsite,
edge fields and the retained escape potential are onsite at edge centers, and
bounded hopping joins incident sites. For a degree-bounded graph family,
interaction-picture finite-range Lieb–Robinson estimates apply to bounded even
local observables, including gauge-invariant densities. The unbounded electric
terms are onsite and can be removed by their local unitary interaction picture;
the norms of the hopping interactions stay bounded. The gauge constraints do not
prevent this proof, since it may be established in the ambient local algebra and
restricted to the invariant physical subspace.

The d>2 crossing-time corollary is a sufficient approximation regime if a supplied
family already has ballistic transport and diameter O(N^(1/d)). It neither
proves transport nor excludes d<=2, and all d>3 and random bounded-degree hosts
remain alternatives. Its interpretation in the synthesis is appropriately
conditional.

### Convention correction found and addressed

For the source pulse exp(-i theta B), the derivative is

    chi = +2 Im <U psi, O U B psi> = -i <[O_H(T),B]>.

An earlier synthesis draft used +i on the commutator while the protected-host
experiment used +2 Im. The root was notified; the synthesis now has the correct
minus sign. Absolute response-error estimates are unchanged. Source convention
must stay consistent in the final integrated manuscript and regenerated output.

## Induced gauge calculation

The Feshbach identity and use of z excluding kappa*Q are correct. It is an
energy-dependent resolvent compression, not automatically an autonomous vacuum
Hamiltonian. For the complete triangle, neutral rewiring vanishes exactly. Every
minimal six-hop path stays outside the vertex vacuum until its final step. The
negative fermion loop product and five negative energy denominators yield the
reported positive sixth-order matrix element. The distinction between a Taylor
coefficient in t and the frozen t=1 matrix element is explicit and necessary.

For the static fixed-graph comparison, filling the edge-derived low band gives
E=Nm-Tr sqrt(m^2+t^2 BB^dagger). The binomial expansion and the oriented-cycle
coefficient in Tr L^ell yield

    J_ell = ell*C_(ell-1)*t^(2ell)/(2m)^(2ell-1) > 0.

The sign is consistent with the exact geodesic controls. The term favors pi
holonomy of an isolated loop in these conventions, not a chosen nonzero electric
flux configuration. Electric flux is conjugate to holonomy. Frustration between
overlapping loops and the absence of a frozen-coupling large-N approximation
must remain explicit. The fourth-order degree-variance penalty is a constructive
organizing term in its controlled limit, with no dimension preference implied.

The omitted-electric-space block comparison, residual resolvent estimate and
min-max/Temple structure are valid. Their floating evaluation is not rigorous
roundoff enclosure, as the derivation states. Full-space positivity of the
eliminated block is supported by a positive margin much larger than ordinary
roundoff, but the tiny displayed residual/Temple widths are not themselves
certified digits. No massless photon, relativistic cone or tensor graviton is
established by a Wilson-loop vertex.

## Integration recommendation

Lead the theory manuscript with the explicit Gauss-constrained Hamiltonian and
one physical chain: matter hopping creates charged flux records; coherent closed
processes create loop amplitudes and gauge interactions; a persistent graph
record supports local matter response with a calculable full-law error. The new
formation kernel and the conditional locality theorem make that chain genuinely
constructive while exposing exactly which link is still missing: appreciable
formation and persistence of the same relational spatial state from the original
neutral preparation.

Present continuum Einstein/Standard-Model terms as matching obligations, not
additional modules already derived from CP-12. A three-dimensional continuum,
a common Lorentzian cone, two tensor helicities, universal stress coupling,
Standard-Model representations and a physical scale calibration remain open.
The present finite predictions are discriminating statements about this model
and potential quantum simulators; they are not yet empirical predictions beyond
established quantum physics.

## Final integrated-manuscript pass: neutral-origin record and coordination

A subsequent review examined `MANUSCRIPT.md`, `protected_geometry/NEUTRAL_RECORD.md`,
`NEUTRAL_RECORD.json`, `neutral_record.py`, and the added static coordination
comparison. This was a source/algebra pass, without rerunning their expensive
evolutions.

The neutral-origin prototype is a valid operational construction with an explicit
intervention. The unforced wait starts from the original neutral ray. The complete
projective record {P_g G, P_g(1-G)} retains every outcome; only its analysis is
conditional on G. The graph-relative receiver is fixed by the recorded initial
ring and is independent of the final graph, avoiding distance reclassification.
The denominator is the initial recorded-event probability. The source exp(-i theta
n_0) commutes with the record, so no final-event quotient derivative is missing.
This establishes a finite recorded-history response, not autonomous record
formation or an unmeasured geometric phase.

The branch-column projection is a contraction into Frobenius norm, so inherited
waiting error does not acquire a square-root-of-branch-count factor. The bilinear
response error correctly uses actual branch-column norms and separate X/Y grouped
errors. Positive-event amplitude intervals, conservative interval division by p_G,
and capping survival intervals at one are sound.

The energy-change estimate properly uses the bounded observable
A=Phi(H)-H=Phi(K+R_off)-(K+R_off). The diagonal unbounded electric contribution
commutes with every record projector. Thus the quoted wavefunction error bounds
the change in energy without attempting to bound an unbounded energy expectation
from norm alone. The complete-record entropy is the Shannon entropy of its
orthogonal branch weights because its input is pure. A finite 30-outcome entropy
continuity bound is appropriate even though the physical rotor space is infinite.
The stated total-variation estimate is conservative. Apparatus preparation and
reset costs remain uncomputed, and global system-plus-memory purity can remain
unchanged. The component note makes these limitations clear.

Two precise wording recommendations were sent to the integrating worker:

1. The prototype's final `Gcirc` projector tests full circulation in either
   orientation. Its record contains graph and event flag, not oriented flux or
   individual field values. Replace the manuscript's 'original circulation
   survival' with 'full circulation probability at final time'. Same-host
   persistence is measured; persistence of the same oriented flux configuration
   is not established.
2. The new coordination test uses all 455 N6,Q12 graphs at **supplied zero
   incidence angles** and m=4,8. Its matrix function is exact, but its evaluations
   and minimizer comparison are floating eigensolves. Describe this as computed
   minima of the exact static-band energy at those fixed angles, rather than
   mathematically certified minima over the full graph-plus-holonomy sector.
   The fourth-order positive coordination penalty remains an analytic result.

Subject to these scope clarifications, the integrated manuscript distinguishes
construction from matching obligations correctly. Its strongest new operational
chain is now neutral unfiltered evolution -> rare charged-ring event -> explicit
costed diagnostic record -> persistent host and matter response. The rarity,
external recording, tiny graph diameter, and absence of size-dependent spatial
formation are stated. This is a useful constructive prototype and should not be
advertised as autonomous spacetime formation.

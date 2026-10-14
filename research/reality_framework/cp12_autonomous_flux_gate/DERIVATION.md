# Autonomous CP-12 flux gate: derivations and precise scope

Frozen law and couplings are unchanged. This packet extends
`4e2b7a876a3a45ed554408ed803b572dd9cab580`, itself an extension of `03256f9`.
External time, canonical quantum mechanics, compact U(1), the supplied Q
sector and the initial ray remain inputs. None of the results below derives
relativistic spacetime, a dimension, a state-selection law or gravity.

## 1. Audited flux tradeoff and attribution

Write M=N(N-1)/2, n_v for occupied vertex fermions, h for holes on present
edge modes, r for removable present edges and k for nonremovable present
edges whose edge fermion remains occupied. The absent edge modes are filled.
Conservation of M total fermions gives h=n_v and

\[
 r+k+n_v=Q.
\]

Gauss law makes the two leg fields sum to zero on occupied edge modes and
to one on hole modes. An occupied, flux-bearing edge therefore has at least
two units of electric square; a hole edge has at least one. With
\(F=(\epsilon_E/2)\sum E^2\),

\[
 F\ge\epsilon_E k+\frac{\epsilon_E}{2}n_v
   =\epsilon_E(Q-r-\tfrac12 n_v),
\]
\[
 2m n_v+F\ge\epsilon_E(Q-r)
            +(2m-\epsilon_E/2)n_v.
\]

All are diagonal operator inequalities on the unbounded physical space.
Rewiring swaps one neutral occupied edge for another and preserves r;
the diagonal energy also preserves r. Thus

\[
 [R,r]=0,\qquad \frac{d\langle r\rangle}{dt}
             =i\langle[H_{\rm hop},r]\rangle.
\]

The independent extension's exact rational checker has now been executed,
not merely read. It passes all 30, 78, 872 and 3975 Gauss states in its four
declared finite samples. Those samples check an algebraic proof; they are
not the justification for extending the inequality to unbounded fields.

## 2. A total-energy constraint that includes negative hopping

Positive component costs alone cannot exclude protection, since hopping
can have negative expectation. The following bound explicitly includes it.
Assume \(a=2m-\epsilon_E/2\ge0\). In the rotor-angle representation, for
each graph the hopping-plus-a-mass quadratic fermion matrix is

\[
 \begin{pmatrix}aI_N&tB(\theta)\\tB(\theta)^\dagger&0_Q\end{pmatrix},
 \qquad \mathop{\rm tr}BB^\dagger=2Q.
\]

There is one unit-modulus entry per half-leg, with the frozen orientation
sign. If s_j are the singular values of B, its negative one-particle
eigenvalues are \((a-\sqrt{a^2+4t^2s_j^2})/2\). Fermionic occupation
therefore bounds the many-body operator from below by their sum. Padding
with zero singular values to N and using concavity of the square root gives

\[
 H_{\rm hop}+a n_v\ge -N C(\nu),\quad
 C(\nu)=\frac{\sqrt{a^2+8t^2\nu}-a}{2},\quad\nu=Q/N.
\]

This pointwise angular bound holds before restricting to Gauss states and
therefore holds after that restriction. It holds on every graph block.
The rewiring Laplacian is positive. Combining the bounds proves

\[
 \boxed{H\ge\kappa Q+\epsilon_E(Q-r)-NC(\nu).}
\]

For the uniform neutral initial ray, the exact conserved energy is
E=\(\kappa Q\), so at every time

\[
 \boxed{\frac{\langle r\rangle}{N}\ge
         \nu-\frac{C(\nu)}{\epsilon_E}.}
\]

At frozen m=t=1 and epsilon_E=1/4, a=15/8. The right hand side is positive
for nu>49/2. Thus low mean removability is impossible at these densities
for this initial energy, even allowing the most favorable negative hopping.
At nu=30 the removable fraction is at least 0.0846675; at nu=50 it is at
least 0.2714920. These are examples of the formula, not large-N evolutions.
The bound is vacuous at nu=1,2,3: it does not exclude the proposed mechanism
in the physically relevant sparse sectors, including a candidate degree-six
geometry. It does not assert sharpness or settle the phase diagram.

## 3. Protection is not a destination preference

Fix the fermion/electric configuration mu. Let P(mu) be the set of potential
edges that are nonneutral (edge hole or nonzero field). All must be present,
by the absent-edge constraint. The remaining present edges are removable.
Within this fixed-mu fiber, R is exactly the Johnson-graph Laplacian on
choosing r edges from K=M-|P| neutral potential edges, with coefficient
Gamma/M. Its degree is r(K-r)=r(M-Q). Its zero eigenvector is the uniform
superposition of the allowed r-subsets. No distance classifier occurs in
this fiber law. For r=0, R vanishes regardless of the geometry of P.

There is also an exact state-dependent destination identity. Suppose an
arbitrary coherent physical state psi has two specified potential edges f,g
absent throughout its support. Define the deletion map D by coherently
removing one removable edge. Directly from the frozen rewiring amplitudes,

\[
 \|n_f H\psi\|^2=\|n_g H\psi\|^2
                 =\frac{\Gamma^2}{M^2}\|D\psi\|^2.
\]

The two target vectors differ only in which of the commonly vacant neutral
bits is activated; their matter/field amplitudes are identical. Hopping
and diagonal terms cannot create an absent edge. The identity includes
coherent interference among all deletion parents.

This concerns instantaneous cross coupling (and its leading small-time
probability), not all finite-time amplitudes. Higher-order matter dynamics
can distinguish topology, as the earlier exact certificates already show.
It must not be promoted to a theorem that no emergent local phase can exist.

An exact prepared comparator has a flux-protected triangle with a neutral
three-edge tail on six vertices. Vacant destinations at initial distances
2,3,4 each have coupling square 1/7500. Protecting the triangle removes
removal opportunities but supplies no near-destination filter in this test.
This comparator is not used as an initial state of the formation runs.

The deletion-map destination identity is inherited from the earlier
cp12_coherent_destinations derivation, rather than claimed as newly found.

If n_v=0, all present edge fermions are occupied and the flux is
divergence-free. Summing Gauss law on either side of a graph bridge proves
that bridge's flux is zero. Vacuum flux protection therefore lives on cyclic
parts of the graph. It does not select long loops, a dimension or locality.

## 4. Protected edges are internal graph records

Use the common kinematic embedding with graph basis g and complete
fermion/electric basis mu. A physical pure state has amplitudes A(g,mu).
Gauss and the absent-edge constraints imply

\[
 A(g,\mu)\ne0\ \Longrightarrow\ P(\mu)\subseteq g,
 \qquad |P(\mu)|=Q-r(\mu).
\]

Thus measuring a complete environment configuration reveals the protected
core, leaving at most

\[
 K(r)=\binom{M-Q+r}{r}
\]

compatible graphs. For the classical joint basis measurement,

\[
 \boxed{H(G\mid\mathrm{configuration})\le
        \langle\log_2 K(r)\rangle.}
\]

This is a support bound, not a claim that the right hand side is attained.
In the r=0 sector the graph is exactly P(mu), so the graph reduced density
matrix is diagonal. There are no cross terms between r=0 and r>0 after
tracing out mu, because r is determined by mu alone in fixed Q.
Complete protection supplies a which-graph record even for a globally pure
state. It does not guarantee that the graph is spatial or that the record
persists under subsequent hopping.

Large graph entropy is also insufficient to infer loss of the dangerous
conditional coherence. The Johnson incidence identity on each fixed-mu
fiber has D^dagger D=r I+A_J and
R=(Gamma/M)[r(K-r)I-A_J]. Since K-r=M-Q,

\[
 \boxed{D^\dagger D=(M-Q+1)r-\frac{M}{\Gamma}R.}
\]

This extends the earlier neutral-pool incidence identity to arbitrary
physical matter/field fibers. It evaluates the actual state's coherent
deletion norm without an inserted graph record. Complete configuration
dephasing would give only mean r; it is a hypothetical, energy-changing
comparator. The amplification is bounded by M-Q+1. The numerator operator
norm is at most Q(M-Q+1), giving its full-law error directly from the state
bound. The saved audit divides bounded numerator/denominator intervals.
This norm is a residual capacity diagnostic, not a matter-response or causal
cone observable. In a state with no universally absent edge it must not be
inserted into a common-vacant formula without the latter's support condition.

The numerical audit unfolds the signed permutation quotient into normalized
physical amplitudes, constructs A as a sparse matrix, and evaluates
rho_G=A A^dagger. Total entropy remains zero. Entropy of rho_G measures
graph/environment entanglement in this stated embedding. The graph bit is
gauge neutral; every expanded state obeys the frozen Gauss constraints.
The r=0 block is inspected as a diagnostic and is not selected or fed back
into the main dynamics. An independent full physical matrix checks the
partial-trace spectrum and all orbit signs.

## 5. Continuous-time bounds, not just sampled minima

On any fixed-Q physical basis state there are at most 2Q hopping transitions,
one per half-leg. Each has magnitude |t| and changes r or k by at most one.
The commutator row and column absolute sums are at most 2Q|t|. The Schur
bound gives

\[
 \|[H,r]\|,\ \|[H,k]\|\le2Q|t|,
 \qquad \left|\frac{d\langle r/Q\rangle}{dt}\right|,
 \left|\frac{d\langle k/Q\rangle}{dt}\right|\le2|t|.
\]

For k, R also preserves the fixed matter/field protected set. Unbounded
diagonal electric energy commutes with these bounded counts, so it does not
invalidate the rate bound. For grid spacing h, every time is at most h/2
from a sampled time. Subtract |t|h and the sampled expectation error from
the minimum to obtain a lower bound over the entire covered interval.
At frozen t=1 and h=0.1, the interpolation allowance is 0.1. We use this
coarse, uniform bound rather than inferring persistence from a plot.

These bounds concern expectations in the unforced state. A positive
probability of an r=0 component is allowed; a large mean r is not a no-go
for a rare conditional branch or for a different physical preparation.
The flux persistence bound concerns its mean amount, not survival of the
same protected edges or loops and not stability of an intrinsic metric.

## 6. Numerical projection, leakage and finite-time state evolution

The retained space has total electric square at most the declared budget.
It is closed under all retained transitions; B keeps every outgoing
full-law transition. The signed permutation quotient is exact. Ordinary
SciPy matrix exponentials evolve the projected Hamiltonian over finite
intervals. For exact U and projected V, Duhamel gives

\[
 \|U(T)\psi-V(T)\psi\|\le\int_0^T\|BV(t)\psi\|dt.
\]

The fifth-order midpoint leakage-vector bound is inherited unchanged from
cp12_distance_history/evolution.py. It bounds a quadrature integral; it is
not a time-Taylor extrapolation for physical observables. Segments' bounds
add. For a normalized positive contraction O and state error e, use the
conservative expectation bound 2e+e^2. Counts r/Q,k/Q,n_v/Q and the
low-removability event are such contractions. Full graph-distribution TV
error is at most e+e^2/2. Norm error alone does not certify the unbounded
electric-energy diagnostic. However, exact conserved total energy supplies
an additional inference. Write H=kappa Q+F+D, where
D=2m n_v+H_hop+R is bounded, with

\[
 \|D\|\le K=2mN+|t|\sqrt{2NQ}
                     +2\Gamma Q(M-Q)/M.
\]

Both the exact and projected evolutions conserve their known mean total
energy kappa Q. Therefore the full-law electric-energy expectation differs
from its projected value by at most K(2e+e^2), plus the saved projected
energy/norm drift terms. ENERGY_AUDIT.json applies this inference; it does
not directly apply a bounded-observable estimate to F. The argument also
guarantees finite electric-energy expectation, since F is total energy
minus a bounded operator.

For the dense-sector count audit, the detector-specific method also applies
to a source-free real expectation. If O is a contraction commuting with
the retained projector, e is the forward error and d is the backward
leakage bound from the projected final vector O x_projected, then

\[
 \boxed{E_O\le2 e d+e^2.}
\]

Expand the expectation in the forward state error. The quadratic term is
at most e^2. Write the linear term using the forward Duhamel residual and
back-propagate O x_projected. The projected backward path is orthogonal to
the outgoing residual, leaving only its error. The cross integral is at
most ed, and the two real cross terms give 2ed. This estimates the unforced
state without inserting a measurement or physical time reversal. Eight
independent dense expectations and a complete physical matrix check it.

For graph entropy, partial trace contracts trace distance and normalized
pure-state trace distance is at most the vector distance e. On graph
dimension G use e log2(G-1)+h2(e) when e<=1-1/G, otherwise log2 G.
Only this graph entropy receives a continuity error in the saved record
audit; the other reduced-state numbers are projected diagnostics.

All bound evaluations and matrix computations use ordinary floating point.
There is no interval-arithmetic certification of accumulated roundoff.
Independent matrices, budget refinement, norm/energy checks and finite
source pulses are numerical cross-checks, not replacements for that caveat.

## 7. Matter response after partial protection

No initial-distance record or geometry postselection is inserted. Evolve
the same uniform neutral ray to declared waits 5 and 10. Apply only the
gauge-invariant diagnostic exp(-i eta n_0). Let O be mean density at all
other vertices, independent of graph-distance bins. Since S=n_0 is a
projector and x=U v, y=U S v,

\[
 f(\eta)-f(-\eta)=4\sin\eta\,\mathrm{Im}\langle x|O|y\rangle,
 \quad \chi=2\mathrm{Im}\langle x|O|y\rangle
          =\frac{f(\pi/2)-f(-\pi/2)}2.
\]

This finite-pulse extraction has no finite-difference bias. Each pulse is
a contraction of the waiting-state error; outgoing residuals bound its
post-source error. Half the sum of the two expectation-error bounds controls
chi. Turning off hopping only after the same full-H waiting preparation is
an explicitly diagnostic comparator: its source commutator is zero.

Pulse work is not free or matched between signs. The source commutes with
R, mass and electric energy. Its work operator is a difference of bounded
hopping operators, with norm at most 2|t|sqrt(2NQ); the same singular-value
argument bounds hopping alone. This makes the saved work-error estimate
legitimate despite unbounded electric fields. The source commutes with r,
so it does not instantaneously select a more protected geometry.

## 8. Initial resources and the actual remaining scientific gate

The initial state is the same neutral equal-amplitude ray over all Q-edge
graphs: energy kappa Q, variance 2Qt^2, zero total entropy, zero flux, r=Q.
Its graph entropy is initially zero even though a graph-basis measurement
has entropy log2 binom(M,Q). It supplies coherence across all these graph
configurations. A short combinatorial description is not a local physical
preparation mechanism. One can describe a sequential counter-controlled
Dicke-state circuit that samples each next graph bit with conditional
amplitude sqrt(q_remaining/K_remaining); its telescoping amplitudes give
the uniform Q-subset state. Such a circuit requires an external controller,
graph registers and resource handling, none derived from the frozen law.
No local gate count, physical cost or autonomous preparation time is claimed.

Main Q=N runs at N=4,5,6 permit a controlled finite-size comparison at common
times, but a connected graph in this sector has cycle rank one. They are
auxiliary mechanism tests, not a thermodynamic 3D gauge-capacity phase.
Additional Q=6,N=5 and Q=12,N=6 controls do not repair that missing large-N
phase test. No spectral/volume dimension is fitted to these small graphs.

Autonomous internal records and persistent partial protection would be real
mechanism-level results. Spacetime emergence additionally requires an
extended spatial family to be selected with controlled size dependence,
stable intrinsic metric records and suppressed distant responses while
retaining local matter propagation. Neither complete protection nor a
finite nonzero susceptibility proves those requirements. The present
frozen law is not modified to force a positive answer.

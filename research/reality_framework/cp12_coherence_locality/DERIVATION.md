# CP-12: coherence, records, and a controlled prepared locality regime

Additive continuation of [`d2fbd05`](https://github.com/ryangrvr/TestingGRUT/commit/d2fbd05c420a209faae618c88fbe4d1902bc999c).
The frozen Hamiltonian and physical Hilbert space are unchanged. Write
\(M=\binom N2\), \(H=H_B+V\), where \(H_B\) retains every graph block's
actual rewiring diagonal and \(V\) contains only graph-changing entries.
The constants are \(m=t=1,\kappa=1/2,\epsilon_E=1/4,\Gamma=1/10\).

The preceding Gauss-law bound and coherent matching trial give, at
\(Q=2N\),

\[
H_B\ge N\ell_N I,\qquad E_0\le U_N,
\]
\[
\ell_\infty=(353-\sqrt{128321})/160,
\quad u_\infty=49/32-\sqrt{21/8},
\quad g_\infty=\ell_\infty-u_\infty
=0.05631513315248\ldots>0.
\]

Here \(U_N\) is the Rayleigh energy of the specified matching tent state.
The inherited lower bound applies to every fixed-graph block by convexity
in the degree, not just regular ones. All analytical statements below
include integer electric fields without a cutoff. Numerical controls
declare their regulators separately.

## 1. Ground-state coherence is required energetically

For any normalized ground state, define

\[
C_{\rm rw}=-\langle V\rangle=\langle H_B\rangle-E_0.
\]

Then \(C_{\rm rw}\ge N\ell_N-U_N\). The exact full operator norm is
\(\|V\|=\Gamma Q(M-Q)/M\), attained in the completely neutral Johnson
fiber. Hence

\[
\boxed{\liminf\frac{C_{\rm rw}}{\|V\|}
\ge\frac{g_\infty}{2\Gamma}=0.2815756657624\ldots .}
\]

This is an invariant energy/coherence statement, not a choice of an
event-labelled order parameter. It remains valid in degenerate ground
spaces and for mixtures supported there.

## 2. Recording geometry has an exact system-energy cost

An ideal, retained orthogonal record of the complete graph implements
\(\mathcal D_G(\rho)=\sum_G P_G\rho P_G\). It preserves the expectation
of \(H_B\) and removes that of \(V\). Therefore

\[
\boxed{\Delta E_{\rm graph\ record}=C_{\rm rw}.}
\]

At \(Q=2N\), \(\liminf\Delta E/N\ge g_\infty\). This is a system energy
increase. An apparatus must exchange energy; no total-energy violation or
autonomous preparation is asserted.

For a weaker record, associate an ancilla with each potential edge:
absent \(\mapsto|0\rangle\), present
\(\mapsto\cos\theta|0\rangle+\sin\theta|1\rangle\). Every rewiring changes
exactly two edge-presence bits, so its record overlap is \(\cos^2\theta\).
Graph-preserving matter/gauge terms have overlap one. Consequently

\[
\boxed{\Delta E(\theta)=\sin^2\theta\,C_{\rm rw}.}
\]

The two single-edge record states have equal-prior trace distinguishability
\(D=|\sin\theta|\). Thus this particular record architecture costs at least
\(D^2g_\infty N+o(N)\) in ground-state system energy. This is not a universal
bound for every measurement architecture or a statement about actual
edge-presence priors.

If \(\Delta E_e\) denotes the energy cost of separately making a perfect
record of edge \(e\), then

\[
\sum_{e=1}^M\Delta E_e=2C_{\rm rw}.
\]

Each graph-changing matrix entry is erased by exactly the two corresponding
single-edge records. The average single-edge cost is \(2C_{\rm rw}/M\),
of order \(1/N\) in the fixed-density ground regime. A complete record
and a single microscopic record have very different resource costs.

`full_ground_records.py` explicitly constructs the system–ancilla pure
state on the four-vertex host and evaluates \(H\otimes I\), independently
checking the overlap result. It also recomputes full ground states and
their graph distributions at five vertices, cutoffs 2 and 3.

## 3. An independently calibrated neighborhood pool

Let \(F\) be a fixed set of \(S\) potential short bonds, specified before
examining the dynamics. Put \(L=M-S\) and let \(k(G)\) count present bonds
outside \(F\). Assume \(Q\le S\). No embedding or spacetime is derived by
introducing this calibration.

In the neutral graph fiber, the permutation quotient of the rewiring
adjacency on levels \(k\) is tridiagonal. Its diagonal is

\[
D_k=(Q-k)(S-Q+k)+k(L-k),
\]

and the off-diagonal is

\[
c_k=\sqrt{(Q-k)(L-k)(k+1)(S-Q+k+1)}.
\]

The two diagonal contributions count short-to-short and long-to-long
swaps. The forward degree between levels is \((Q-k)(L-k)\); the reverse
degree is \((k+1)(S-Q+k+1)\). The square root follows from normalization
of uniform layer states.

For \(P_K=1_{k\le K}\), the neutral internal adjacency norm is exactly
the largest eigenvalue of this quotient truncated at \(K\). Perron
positivity and the short/long permutation symmetry justify the reduction.
All other physical rewiring fibers are submatrices of the neutral graph
matrix, and the neutral fiber is itself physical. Thus these norms are
also the complete physical operator norms. In particular,

\[
\boxed{\|P_KH(I-P_K)\|^2
=\frac{\Gamma^2}{M^2}(Q-K)(L-K)(K+1)(S-Q+K+1).}
\]

Only the rewiring operator crosses this projection. This exact result
includes arbitrary electric fields and matter occupations.

## 4. A fixed sparse neighborhood cannot dominate the ground state

Take \(Q=2N\), \(S=sN+o(N)\) with finite \(s\ge2\). For \(K=\delta N\),
\(0\le\delta<2\), the internal off-diagonal norm obeys

\[
\|P_KVP_K\|\le\Gamma K+O(\sqrt N)+O(1).
\]

The long-to-long diagonal is at most \(\Gamma k\), the short-to-short part
is \(O(1)\), and the cross-layer quotient norm is \(O(\sqrt N)\). The
cross-boundary norm squared is \(O(N)\) at fixed positive \(\delta\).

If \(\Gamma\delta<g_\infty\), the energy in this sector is at least
\(N\ell_N-\Gamma\delta N-O(\sqrt N)\), separated extensively from the
full ground upper bound. The projected ground-state equation yields

\[
\boxed{\Pr_{\rm ground}[k\le\delta N]=O(N^{-1})\to0
\quad\text{if }\delta<g_\infty/\Gamma.}
\]

For \(K=o(N)\), the sharper estimate is
\(\Pr[k\le K]=O((K+1)/N^2)\). Exact quotient/Schur bounds give the finite
certificates in `locality_pool_results.json`; the certificates use outward
Decimal bounds, not a large-Hamiltonian eigensolve.

A mean-value consequence avoids a sector cutoff. Split \(V\) into
short-preserving, long-preserving, and crossing operators. For any state,

\[
-\langle V\rangle\le\Gamma\langle k\rangle+
\frac{\Gamma S^2}{4M}+\frac{2\Gamma}{M}\max_k c_k.
\]

Combining with section 1 gives

\[
\boxed{\liminf\frac{\langle k\rangle}{Q}
\ge\frac{g_\infty}{2\Gamma}=0.2815756657624\ldots .}
\]

It concerns **any independently specified** \(O(N)\)-size neighborhood
pool. It does not exclude a graph-dependent relational metric, or a
reference sector interacting with the geometry and changing the total
Hamiltonian. Unioning over all label permutations or defining the pool
after a graph observation is not this theorem. A coordinate-free spatial
phase is not ruled out by this fixed-calibration statement.

As a diagnostic corollary, an external fixed-pool bias
\(H_\eta=H-\eta N_F\) cannot concentrate its ground state in \(k=o(N)\)
if \(\eta<g_\infty/2\approx0.02815757\). The symmetric matching trial has
\(\langle N_F\rangle=QS/M=O(1)\), while the local sector gains at most
\(\eta Q\). The extensive gap remains. This is only a necessary threshold
for that specified bias; crossing it proves neither localization nor
spatial selection. No bias is added to the frozen model or used in the
executed dynamics.

## 5. Prepared locality has a different time scale

Now prepare any physical state entirely in \(P_0\), with every connection
in the independently specified pool. Matter and gauge Hamiltonian terms
preserve \(k\). Rewiring is the only term changing it. The raising block
from \(k\) to \(k+1\) has norm \(\Gamma c_k/M\), which is at most

\[
g\sqrt{k+1},\qquad g^2=\Gamma^2QLS/M^2.
\]

Writing \(p_k=\|\Pi_k\psi\|^2\), Cauchy-Schwarz gives

\[
|d\langle k\rangle/dT|\le
2g\sum_k\sqrt{k+1}\sqrt{p_kp_{k+1}}
\le2g\sqrt{\langle k\rangle}.
\]

Starting at zero and regularizing the square root if necessary gives

\[
\boxed{\langle k(T)\rangle\le g^2T^2.}
\]

For \(Q=\nu N\), \(S=sN\), \(g^2\to2\Gamma^2\nu s\). Markov's
inequality implies \(\Pr[k\ge\varepsilon N]\le g^2T^2/(\varepsilon N)\).
A macroscopic density of outside bonds therefore cannot develop on
\(T=o(\sqrt N)\), even for arbitrary coherent preparations and full
matter/gauge dynamics. This is a topology-density bound; by itself it
does not prevent a small number of harmful shortcuts or prove a causal
response cone.

## 6. The full detector response approaches local dynamics

There is a stronger controlled regime. Let \(S=Q+r_N\), with a
subextensive number \(r_N\) of pool vacancies. The exact leakage norm from
\(P_0\) is

\[
b_0^2=\frac{\Gamma^2}{M^2}Q(M-S)(r_N+1).
\]

The internal graph-changing norm is exactly
\(v_{\rm in}=\Gamma Qr_N/M\). Let \(A,B\) be bounded physical observables
commuting with graph projections, with \(\|A\|\le1\) and
\(\|P_0BP_0\|\le1\), and let the initial state be any density operator
supported in \(P_0\). Compare full evolution to
the complete independent graph blocks, retaining their actual rewiring
diagonal. Then

\[
\boxed{|R_{\rm full}(T)-R_{\rm blocks}(T)|
\le4b_0^2T^2+4v_{\rm in}|T|.}
\]

Proof: the compressed full propagator differs from the propagator of
\(P_0HP_0\) by at most \(b_0^2T^2/2\), and its complementary amplitude
is at most \(b_0|T|\). A block-diagonal detector's compressed Heisenberg
error is at most \(2b_0^2T^2\); its commutator with the source adds a
factor two. Duhamel comparison of \(P_0HP_0\) to independent graph blocks
adds at most \(4v_{\rm in}|T|\) to the response. No incoherent graph
mixture or initial energy eigenstate is required.

Only the source's prepared-sector norm is needed: because \([B,P_0]=0\),
\(P_0[A(T),B]P_0=[P_0A(T)P_0,P_0BP_0]\). Its larger norm on remote
graph sectors does not enter this response estimate. This lets an
apparatus modulate the full frozen hopping term at one vertex, without
querying which potential bonds belong to the calibrated pool. The
presence factor in the original Hamiltonian activates the existing
bonds. On the prepared sector its norm is at most the square root of
that vertex's pool degree; only a scalar strength calibration is used.

At fixed \(\nu\), the error is

\[
O\big((r_N+1)T^2/N+r_N|T|/N\big).
\]

Thus for \(r_N=o(N)\), fixed-time local detector statistics converge to
the block-local theory despite coherent geometry fluctuations. This
statement is stronger than protection of one initially definite graph.

There is also a back-reaction limitation in precisely this window. For
any bounded graph-only readout \(A_G\), independent graph-block evolution
commutes with \(A_G\), and a graph-diagonal matter source \(B\) commutes
with it too. Thus \(R_{\rm blocks}^{A_G,B}=0\) exactly and

\[
\boxed{|R_{\rm full}^{A_G,B}(T)|
\le4b_0^2T^2+4v_{\rm in}|T|.}
\]

Whenever the matter-response error vanishes by the stated scaling, this
linear graph-connectivity response also vanishes. The controlled regime
therefore supplies propagation on a stable background, not a finite
gravitational response of that background to a bounded matter source.
Different observation times, extensive sources, or a metric encoded in
other physical degrees of freedom need their own calculation. The
connection-presence readout is tested directly in `prepared_response.py`.

## 7. Conditional commutator cone and its physical assumptions

Assume the pool's incidence graph has maximum degree \(D\) uniformly
bounded. Group each edge's presence bit, edge fermion, and two compact
rotors into one edge-center cell; each vertex has its fermion cell.
In an independent graph block:

* mass is on-site at vertex cells;
* electric energy and the removable-edge diagonal are on-site at edge cells;
* every hopping interaction joins an incident vertex and edge cell and
  has norm at most \(|t|\).

Removing absent bonds reduces the interaction set. Each term is even
in fermions. Gauge-invariant sources and detectors preserve the physical
Gauss sector. Disjoint even observables commute; a spin-chain
Jordan-Wigner string must not be mistaken for physical nonlocality.

Put the unbounded on-site electric Hamiltonian into the interaction
picture. It preserves cell supports and hopping norms. The usual
interaction-chain recursion now applies to bounded, time-dependent even
interactions. Primary antecedents are Nachtergaele and Sims,
[arXiv:1410.8174](https://arxiv.org/abs/1410.8174), Theorem 3.1, for
unbounded on-site terms; and Nachtergaele, Sims and Young,
[arXiv:1705.08553](https://arxiv.org/abs/1705.08553), Theorem 3.1, for
fermionic locality. We apply their support-growth argument here, not a
new relativistic principle.

For completeness, a conservative direct chain estimate on the incidence
graph is

\[
\|[A_X(T),B_Y]\|\le
\frac{2D|X|}{D+1}\|A\|\|B\|
\sum_{n\ge d_F(X,Y)}\frac{[2|t|(D+1)|T|]^n}{n!}.
\]

At most \(D|X|\) hopping terms meet the initial support. Each subsequent
term has at most \(D+1\) overlapping hopping terms: the incidence edge's
vertex has at most \(D\) legs and its edge-center has two. Including
the term itself overcounts safely. A chain reaching \(Y\) has length
at least \(d_F(X,Y)\). Iteration supplies \((2|t||T|)^n/n!\). Even CAR
operators plus rotor factors retain the same disjoint-support argument.
Uniform finite-cutoff estimates pass to the on-site strongly continuous
dynamics; bounded physical probes remain bounded on the Gauss sector.

The tail also obeys, for any \(\mu>0\),
\(\sum_{n\ge d}x^n/n!\le\exp(xe^\mu-\mu d)\). Combining with section 6
gives an exponential cone term **plus** the explicit full-model error.

This is a Lieb-Robinson-type approximate commutator cone in the supplied
incidence metric. It is not a strict zero outside a Lorentzian light cone,
does not derive an invariant speed, and also holds on nonspatial
bounded-degree hosts. The external time in CP-12 remains an input.

`light_cone_control.py` tests operator norms on complete conditional
path/cycle blocks. Independently, integer nested commutators vanish
exactly through every order shorter than the declared interaction chain.
`prepared_response.py` tests the **full** Hamiltonian with a local frozen
hopping source and vertex-density detectors, including independent
finite-impulse derivatives. No geometry-role oracle enters these probes.

## 8. A macroscopic three-dimensional prepared benchmark

Choose a supplied periodic diamond-type graph with two vertices per
primitive cell and \(N=2L^3\). Each A vertex connects to B vertices in
its own cell and its three preceding coordinate cells. It is connected,
four-regular, and has \(Q=2N\). Optional \(r_N\) extra disjoint short
A–A bonds form a degree-at-most-five pool of size \(S=Q+r_N\); the
physical fixed-Q preparation leaves \(r_N\) vacancies coherently.

The supplied graph's Bloch Laplacian bands are
\(4\pm|1+e^{ik_x}+e^{ik_y}+e^{ik_z}|\). The lower band is a positive
three-dimensional quadratic form near zero. Sparse graph eigenvalues
independently agree with this expression. This checks the **benchmark**
dimension; it is not a CP-12 dimensional-selection calculation or its
complete matter dispersion relation.

At observation times \(T\sim L\), section 6's error tends to zero if
\(r_N=o(L)\). For example \(r_N=\lfloor\sqrt L\rfloor\) permits an
increasing number of coherent geometric vacancies and macroscopic
propagation times simultaneously. In a supplied \(d\)-dimensional host,
\(T\sim N^{1/d}\), the sufficient condition is
\((r_N+1)N^{2/d-1}\to0\). It is compatible with \(d>2\); failure of this
sufficient bound at \(d\le2\) does not exclude those dimensions or select
three dimensions.

One physical preparation is the uniform superposition of neutral
\(Q\)-bond graphs inside the pool, with vertex fermions empty, edge
fermions filled and all fields zero. It satisfies Gauss law exactly and
requires no tailored energy spectrum. Its full-law moments are

\[
\langle H\rangle=\kappa Q+\frac\Gamma M Q(M-S),\qquad
\operatorname{Var}(H)=2t^2Q+b_0^2.
\]

Thus energy per vertex tends to 1.2 and relative energy width tends to
zero like \(N^{-1/2}\). This is a finite-energy-density excited preparation,
not the ground phase. Its coherent pool and physical preparation mechanism
are still supplied. The small full-law matrix control independently
verifies these moments.

## 9. What has and has not been earned

There is now a controlled **prepared** regime in which coherent graph
fluctuations, physical Gauss constraints, and local source–detector
response coexist, extending to macroscopic observation times in a
supplied 3D host. There is also a quantitative explanation for why the
same model's ground state resists a fixed sparse calibration, and for
the dynamical cost of retained geometry records.

The actual formation and selection of the host have not been derived.
Connectedness improves in the full five-vertex ground control, but this
does not prove polynomial spatial growth. A coupled dynamical reference
system is not provided by a fixed pool or an ideal record isometry.
Vertex addressing and scalar source calibration are supplied probe
specifications. If event relabelling is taken as gauge, physically
localized probes additionally require relational matter tags/dressing;
the labelled calculation does not supply them.
Neither Lorentz invariance, a gravitational constraint, universal spin-2
coupling, observed particle physics, nor a complete theory of reality
follows from this packet. The next construction must confront those
obligations rather than reinterpret a conditional propagation theorem
as their solution.

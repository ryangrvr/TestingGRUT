# CP12: a conditional operational cone on prepared bounded-degree geometry

Author derivation, pending independent review. Frozen LAW.md at 0014d615 is unchanged. hbar=1 and external time T are supplied.

## Exact comparator, including the rewiring diagonal

Fix a graph G with Q present connections and M=N(N-1)/2 possible ones. Let P=P_G include every permitted matter/electric state on that graph. Let H_P=PHP restricted to P. It is NOT the Hamiltonian obtained by deleting the entire graph Laplacian.

Let R_e be the projector that edge e is filled and both electric half-fields are zero. The degree of a state in the reconnection graph is (M-Q) sum_(e in G) R_e. Therefore

 H_P = H_bare,G + a_N sum_e R_e
     = H_bare,G + a_N Q I - a_N sum_e (I-R_e),
 a_N = Gamma(M-Q)/M.

The common a_N Q may be subtracted; the last term is a bounded, edge-local matter/electric potential. Its coefficient tends to Gamma, not zero, when Q proportional to N. Thus even a frozen graph retains a finite dynamical consequence of the positive graph Laplacian. The partner audit correctly defines this compressed comparator in its complete report.

## Quadratic detector bound

Write H in P and P-perpendicular blocks, with C=P-perpendicular H P. Different single swaps have orthogonal final geometries and each has norm at most Gamma/M, hence

 ||C|| <= b_N = (Gamma/M) sqrt[Q(M-Q)].

For any initial vector psi in P with norm<=1, define x(T)=P exp(-iHT)psi, y(T)=P-perpendicular exp(-iHT)psi, and x_0(T)=exp(-iH_P T)psi. Variation of constants gives

 ||y(T)|| <= b_N |T|,
 ||x(T)-x_0(T)|| <= b_N² T²/2.

The second bound follows by integrating C†y, not by discarding virtual return paths. It holds for any self-adjoint diagonal blocks because only the off-block perturbation is bounded; the unbounded electric operator does not invalidate it.

If O is bounded, ||O||<=1, and [O,P]=0, its expectation has no x-y interference term. Consequently

 |<O>_full(T)-<O>_P(T)| <= 2 b_N² T².

Indeed the difference of x-sector quadratic forms is at most 2||x-x_0||, and the y-sector contribution is at most ||y||². Equivalently the compressed Heisenberg-operator difference has norm<=2 b_N² T². The estimate includes bounded graph-dependent detectors that preserve geometry; it is not valid for a detector that coherently changes geometry.

For bounded Hermitian source B with [B,P]=0 and ||B||<=1, define the retarded expectation R_AB(T)=i<psi|[A(T),B]|psi>. If A also preserves P and has norm<=1,

 **|R_full(T)-R_P(T)| <= 4 b_N² T².**

No stationarity of psi, electric cutoff, small-T expansion, Gaussian approximation or energy eigenstate is required. Mixed states supported on P satisfy the same bounds. The inequalities become vacuous at sufficiently long times; they are not convergence predictions at T proportional to sqrt(N).

For Q=nu N, b_N²~2 Gamma² nu/N. The observable correction is O(T²/N), sharpening the audit's valid O(T/sqrt(N)) estimate for its more general detector class.

## Local propagation of the compressed theory

Use the bipartite incidence graph of G: each original vertex is one fermion site; each present edge center is one site holding its edge fermion and two electric rotors. The mass, electric energy and a_N R_e are ON-SITE terms. Electric energy is unbounded but self-adjoint. Fermion hopping connects adjacent vertex/edge-center sites, is even under fermion parity, and has bounded norm <=t_hop because rotor U is unitary. Cutoff-projected U has norm<=1 as well.

For a family of prepared graphs with uniformly bounded degree, these are bounded finite-range interactions of uniformly bounded multiplicity. Passing to the interaction picture of the on-site terms preserves support, parity and hopping norm. The usual commutator expansion for parity-even observables therefore gives Lieb–Robinson bounds with constants independent of N and electric cutoff. One may work first in the unconstrained active-site Hilbert space and then restrict to the invariant Gauss/number sector; gauge-invariant observables and dynamics preserve it, so projection cannot increase the bound.

For bounded gauge-invariant EVEN local A and B, with uniformly bounded supports X,Y at incidence distance d_G(X,Y), there are finite constants C_0,a,v_LR independent of N such that

 |R_P(T)| <= C_0 exp[-a(d_G(X,Y)-v_LR |T|)].

The same follows for finite local interventions from the norm commutator estimate. Neither the lattice graph nor a Lieb–Robinson velocity is a Lorentzian metric. The inequality is an upper bound on information propagation, not a dispersion relation or proof of ballistic transport.

Combining with the exact quadratic comparison establishes

 **|R_full(T)| <= C_0 exp[-a(d_G-v_LR |T|)]
                  + 4 Gamma² T² Q(M-Q)/M².**

This is a prepared-state operational response bound for the full frozen graph-changing Hamiltonian, not a global operator-norm cone on arbitrary fluctuating geometries. It controls virtual shortcuts as well as final graph leakage. Source and detector operations are mathematically local, gauge invariant and geometry preserving; a complete autonomous reference apparatus is still not constructed.

## Growing-time window and dimensional qualification

If T_N=o(sqrt(N)) and d_N-v_LR T_N tends to infinity, the right-hand side tends to zero. Thus the model has a controlled growing-time, growing-distance conditional regime; microscopic shortcuts do not invalidate every effective-locality statement.

For a SUPPLIED approximately D-dimensional bounded-degree family, linear size scales as N^(1/D). At times T_N=c N^(1/D), the graph-change correction scales as N^(2/D-1): it vanishes for D>2, is only O(1) at D=2, and this bound gives no useful whole-system control at D=1. These are times proportional to diameter, not a proof that any actual excitation crosses the sample ballistically. Dimension is supplied, not selected; all D>2 share the favorable exponent. This is NOT an argument uniquely selecting three dimensions.

There is a real tension: the previous canonical entropy theorem can make these very spatial families atypical at fixed positive temperature. Conditional preservation of a prepared geometry is different from forming it or making it thermodynamically typical.

## Hostile cases and scope

Unbounded-degree hubs have no N-independent velocity guarantee from this argument. A growing number of initially coherent graph sectors is not covered by the quadratic theorem. In fact two neighboring graphs with a relative imaginary amplitude carry an O(T) probability current, checked explicitly in operator_controls.py. For superpositions supported on a fixed K graphs, the more conservative all-block-diagonal comparison satisfies a state-distance bound <=sqrt(K)b_N|T| by the triangle inequality; arbitrary broad coherence defeats that bound's usefulness.

At T=c sqrt(N) the error bound is O(1); no geometry-formation rate or transport separation follows. At T of order sqrt(N)/Gamma the departure allowance is also order unity. Numerical tests at T=sqrt(N) or 2sqrt(N) with Gamma=.1 are well below that latter scale. Nothing here establishes Lorentz invariance, a graviton, gravitational universality or the cosmological preparation of G.

## Literature and mathematical status

The fixed-geometry ingredient is established theory, not a GRUT novelty: Nachtergaele and Sims, *On the dynamics of lattice systems with unbounded on-site terms in the Hamiltonian*, arXiv:1410.8174 (2014), supplies the relevant unbounded-on-site framework. Nachtergaele, Sims and Young, *Lieb–Robinson bounds, the spectral flow, and stability of the spectral gap for lattice fermion systems*, arXiv:1705.08553 (2017/2018), supplies the parity-even fermionic version. Their ingredients apply here through the explicit interaction-picture/support argument above. The model-specific combination with the CP12 off-block bound is an author deduction; independent proof review and priority assessment remain open.

## Extension: classical geometric uncertainty is allowed

Let D(rho)=sum_G P_G rho P_G denote dephasing between graph sectors and H_diag=sum_G P_G H P_G. For ANY initial state satisfying rho=D(rho), convexity extends the quadratic bound to the entire mixture:

 |R_full(rho,T)-R_diag(rho,T)| <= 4 b_N² T².

There is no restriction on how many graph sectors the incoherent mixture occupies, provided Q is fixed and source/detector preserve each sector. This does not follow for a coherent superposition, but it does cover classical uncertainty about geometry.

If the mixture assigns probability 1-p_bad to bounded-degree geometries where the calibrated source/detector supports have distance at least d, then

 |R_full(rho,T)| <= C_0 exp[-a(d-v_LR |T|)] + 2 p_bad + 4 b_N² T².

The 2 p_bad term is the trivial commutator bound on the remaining geometries. For a general state define delta_coh=||rho-D(rho)||_1. Since the response operator has norm<=2, adding **2 delta_coh** gives a valid bound for the general state's full response. This is a sufficient condition, not a necessary criterion for all coherent states; some large-coherence states can still be local.

The result separates three physically distinct requirements: the weight of unfavorable geometries, the initial graph coherence, and the rate of topology changes. Decoherence may reduce the second; it does not automatically improve the first. No apparatus or cosmological process preparing such a mixture has been supplied by this argument. mixture_control.py explicitly checks quadratic changes for a two-graph mixture against the first-order current in its coherent counterpart.

### The incoherent preparation is not an energy-free operation

For the equal-amplitude neutral filled-edge superposition over all Q-edge graphs, the reconnection Laplacian annihilates the state. Completely dephasing its graph label raises the CP12 system's mean energy by exactly

 Delta E_dephase = Gamma Q(M-Q)/M ~ Gamma nu N.

The fermion, mass and electric diagonal expectations are unchanged. A single neutral definite graph has the same reconnection expectation as the incoherent mixture. This is an extensive but finite energy cost per vertex, not an impossibility theorem or a claim about the apparatus's total energy. It prevents silently substituting an incoherent prepared geometry for a coherent low-energy state at no cost. The N=4,Q=4 value 2/15 is checked in mixture_control.py. Preparation, energy exchange and favorable geometry weights remain separate physical obligations.

# CP-12 relational record gate: derivation

14 October 2026. Additive continuation of `a351fea238f819f29ea1edbd8b7559d43b0be732`, whose parent is `8ece3e20a331b8fe0c745d6ab1d778185e56dade`. The frozen law and supplied neutral initial ray are unchanged.

## 1. Protected cores are electric-support records

Write a physical complete configuration as \(s=(g,b,E)\). There are \(M=\binom N2\) potential edges and exactly \(Q\) present graph edges. An absent edge has its edge fermion occupied and both fields zero. A present edge is removable precisely when its fermion is occupied and both fields vanish.

The edge Gauss equation is

\[
E_{ie}+E_{je}=1-n_{b_e}.
\]

A hole therefore forces nonzero electric field. On the physical Hilbert space the protected-edge projector can consequently be written

\[
p_e=1-\Pi_{E_{ie}=0,E_{je}=0}.
\]

Its support is necessarily a present graph edge. The protected core \(P(s)\) is the unlabelled graph formed by these edges, retaining all \(N\) vertices, including isolates. Its size is \(p=Q-r=k+n_v\). Here \(k\) counts occupied-edge-fermion edges carrying flux, while \(n_v\) counts vertex fermions and, equivalently, present edge-fermion holes. The core is a gauge-invariant final-configuration diagnostic. This definition does not introduce a measurement during evolution or select a classical outcome.

Neutral rewiring only swaps a removable present edge with an absent edge, keeping \(b,E\) fixed. Hence, as full physical operator identities,

\[
[H_R,p_e]=[H_D,p_e]=0,
\qquad
\boxed{[H_R,F(P)]=[H_D,F(P)]=0}
\]

for every bounded function of the protected core. In particular,

\[
\boxed{\frac{d}{dT}\langle F(P)\rangle
=i\langle[H_{\rm hop},F(P)]\rangle.}
\]

This extends conservation of \(r\) under rewiring to the complete protected-edge configuration and its topology. It does not remove the influence of rewiring on later dynamics: rewiring changes the state on which hopping subsequently acts. It also does not assert that full-graph observables commute with rewiring.

## 2. Exchangeability

The physical vertex-permutation unitaries include fermionic signs and the electric-leg relabelling. The initial density operator is invariant and the frozen Hamiltonian commutes with those unitaries. Edge transitivity and the exact sums \(\sum_e n_e^G=Q\), \(\sum_e p_e=Q-r\) imply

\[
\langle n_e^G\rangle=Q/M,\qquad
\langle p_e\rangle=(Q-\langle r\rangle)/M.
\]

The quotient used here is the existing exact signed permutation sector, not a declaration that permutations have been derived as gauge redundancy. Any diagonal label-invariant statistic is constant on its physical orbit. If a normalized orbit amplitude is \(v_j\), its total configuration probability is \(|v_j|^2\), regardless of the orbit size \(s_j\). Each labelled configuration in that orbit has probability \(|v_j|^2/s_j\).

## 3. Geometry and reference distributions

The catalogue enumerates every simple graph on \(N\le6\), partitions it into vertex-permutation orbits and computes component sizes, cycle rank, cycle vertices, triangles, squares, short-cycle packing and intrinsic balls. All reported scores lie in \([0,1]\). Ball volumes average over every root, including isolates; unreachable vertices never enter a finite-radius ball.

Four conditional references are used.

1. **Uniform core size:** at the observed \(p\), distribute the core uniformly among all \(\binom Mp\) labelled graphs.
2. **Parent thinning:** at the observed full graph \(g\) and core size \(p\), choose \(p\) edges uniformly among the \(Q\) edges of \(g\).
3. **Gauss/resources:** among physical configurations represented by the existing reachable signed-invariant basis, choose uniformly conditional on \((n_v,E^2,r)\), with \(E^2=\sum_\ell E_\ell^2\).
4. **Gauss/parent/resources:** additionally condition on the unlabelled parent graph. Representatives have canonical graph masks, so that mask can be used as a class key.

The first two are geometric distributions; arbitrary protection assignments need not admit the stated Gauss-law resources. The last two use actual physical configurations but only the represented reachable orbit support, excluding sign-forbidden or unrepresented configurations. They are not asserted to be the full physical microcanonical ensemble. A configuration-diagonal mixture over those orbits is permutation invariant as a density operator but need not retain the initial pure ray's representation-sector support or its purity.

For either Gauss reference let \(a(j)\) be its conditioning group. Its orbit weights are

\[
K_{ij}=\mathbf1_{a(i)=a(j)}\frac{s_i}{\sum_{\ell:a(\ell)=a(j)}s_\ell}.
\]

This preserves the complete group-probability distribution, not just its mean. Thus \(n_v,E^2,r,k,p\) and the diagonal mass, electric and physical rewiring-diagonal energy components are matched. The parent-conditioned version also matches the full-graph shape distribution.

For application to the full unbounded-electric-field evolution, define this kernel to be the identity on invariant configuration orbits outside the declared represented retained set. This is an explicit bounded extension, not a claim to have enumerated infinitely many configurations. Contributions from those outside orbits cancel in the actual-minus-reference shape distribution. The Gauss-reference results therefore concern this specified finite-support randomization, with its explicit extension; they must not be generalized to every possible resource-matched null.

If \(f\) is a score, its reference expectation is \(\sum_j |v_j|^2(K^Tf)_j\). Reference scores are contractions. The contrast is the expectation of the single diagonal observable \(f-K^Tf\), whose norm is at most one. It is not necessary to add two independent error bars.

For core isomorphism classes \(c\), let \(a_c\) and \(b_c\) denote actual and reference probabilities. We report

\[
D_{\rm TV}=\tfrac12\sum_c|a_c-b_c|.
\]

The class labels are canonical masks used as reproducible identifiers; they do not assign physical meaning to vertex labels.

## 4. A sharp obstruction to a classical energy-matched null

Set \(a=\Gamma(M-Q)/M\). The physical diagonal matrix element of the frozen Hamiltonian is

\[
H_{ss}=\kappa Q+2m n_v+F+a r,
\qquad
F=\frac{\epsilon_E}{2}E^2.
\]

There are exactly \(r(M-Q)\) neutral rewiring neighbours. Hopping has no complete-configuration diagonal element. The flux tradeoff gives

\[
F\ge\epsilon_E(Q-r-n_v/2),\qquad r\le Q-n_v.
\]

For the frozen couplings \(\Gamma=0.1,\epsilon_E=0.25,m=1\), we have \(a\le\epsilon_E\). Therefore

\[
\begin{aligned}
H_{ss}-\kappa Q
&\ge\epsilon_EQ+(a-\epsilon_E)r+(2m-\epsilon_E/2)n_v\\
&\ge aQ+(2m+\epsilon_E/2-a)n_v\\
&\ge aQ.
\end{aligned}
\]

Consequently every density operator diagonal in the complete graph/fermion/electric configuration basis obeys the sharp bound

\[
\boxed{\langle H\rangle_{\rm diagonal}
\ge\kappa Q+\frac{\Gamma Q(M-Q)}M.}
\]

Equality is attained by any neutral basis configuration. If \(0<Q<M\), this exceeds the initial pure neutral sea's energy \(\kappa Q\). Thus a fully configuration-incoherent energy-matched null does not exist at that mean energy. The claim concerns complete configuration diagonality, not diagonality of the graph reduced state or decoherence of one subsystem. At \(Q=M\), the strict gap closes.

The displayed Gauss-resource distributions are diagnostics of configuration probabilities. Their diagonal physical mixtures have explicitly reported higher energy. Any coherent state with those same probabilities would require a separate phase and energy construction; none is assumed here.

For an actual mean-energy-matched physical comparator, diagonalize the existing \(N=4,Q=4\) retained Hamiltonian \(H_P\). With its spectral projectors \(\Pi_\alpha\), use

\[
\rho_{\rm ref}=\sum_\alpha
\Pi_\alpha|\psi_0\rangle\langle\psi_0|\Pi_\alpha.
\]

This positive density operator is supported on physical Gauss states. It preserves all projected energy probabilities. Moreover,

\[
\operatorname{Tr}(H\rho_{\rm ref})
=\operatorname{Tr}(H_P\rho_{\rm ref})
=\langle\psi_0|H_P|\psi_0\rangle=\kappa Q.
\]

It is generally coherent in the configuration basis and mixed globally. It is not claimed to match instantaneous resources, total purity, all full-law energy statistics, or full-law stationarity. For example the full energy second moment includes the additional positive term \(\|BX\|_F^2\), where the columns of \(X\) are the dephased block vectors. The comparator is a mathematical control, not a derived physical preparation procedure.

The same eigendecomposition also supplies a maximum-entropy reference within the declared retained signed-invariant sector, \(\rho_\beta=e^{-\beta H_P}/Z\). Choose \(\beta\) so that \(\operatorname{Tr}(H_P\rho_\beta)=\kappa Q\). This positive physical state's full-H mean is again exactly \(\kappa Q\). Its entropy is maximal under that mean-energy constraint **within this finite retained sector**, not within the full unbounded physical Hilbert space. It need not match energy variance, resources or purity. No thermal preparation or new dynamics is asserted.

## 5. Intrinsic candidate classes

For each configuration, distances are computed in \(P\). Edges outside \(P\) fall into three disjoint bins: finite core distance two, finite core distance at least three, and endpoints in distinct core components. The last bin is not assigned an arbitrary large finite distance.

Let \(K_c(P)\) be the number of possible noncore edges in class \(c\), and \(m_c(g,P)\) the number of present removable edges in it. Uniform completion of that same core by \(r\) neutral edges gives

\[
\mathbb E_{\rm completion}m_c
=r\frac{K_c(P)}{M-|P|}.
\]

We compare \(m_c/Q\) with this value divided by \(Q\), defining both to vanish when no noncore edge is available. All actual scores, reference scores and their contrasts are bounded contractions. Candidate pools are selected from the final configuration; they are not fixed initially empty pools, and the earlier common-vacant-edge theorem is not applied to them.

`shortcut_controls.py` also applies each Gauss/resource kernel directly to the joint-configuration score \(m_c(g,P)/Q\). This additionally controls parent topology in the parent-conditioned version. It is distinct from uniform completion. Its bounded contrast and explicit identity extension obey the same error inequality as the core-shape contrasts.

The bare neutral rewiring matrix element is still \(\Gamma/M\) for every allowed destination. Protection removes destinations without adding a core-distance factor. Instantaneous occupation contrasts can depend on coherent history; they do not prove that bare rewiring itself favors one bin. Since completion expectations are functions of \(P\), they commute with rewiring. Thus the rewiring derivative of a contrast equals that of \(m_c/Q\).

For each physical rewiring transition, \(|\Delta(m_c/Q)|\le1/Q\). Its commutator has absolute row sum at most

\[
\frac\Gamma M\frac{r(M-Q)}Q\le a.
\]

Symmetry and the Schur bound give \(\|i[H_R,m_c/Q]\|\le a\). Core-shape derivative norms are bounded by \(2\|H_{\rm hop}\|\le2|t|\sqrt{2NQ}\).

## 6. Full-law errors and the dense TV witness

Let \(v=e^{-iH_PT}\psi_0\), embedded in the full physical Hilbert space, and \(\|e^{-iHT}\psi_0-v\|\le e\). The outgoing matrix \(B=(1-P)HP\) and the existing analytic leakage sums supply \(e\). This audit reconstructs the identical projected state with sparse exponentials and verifies published counts at every already computed time. The mathematical leakage bound depends on the exact projected state, not on which sparse exponential call reconstructs it.

For a Hermitian contraction the expectation error is at most

\[
\delta=2e+e^2.
\]

The configuration-probability TV error is at most \(e+e^2/2\). Applying the same Markov reference to both distributions then bounds the error in their core-shape TV by \(\delta\). Core/current derivative errors use the operator norms above multiplied by \(\delta\).

For a Hermitian contraction \(O\) commuting with \(P\), the earlier source-free dual estimate gives

\[
|\langle O\rangle_H-\langle O\rangle_P|
\le2e\,d+e^2,
\]

where \(d\) bounds backward propagation error from \(Ov\). The cross term uses orthogonality of a retained backward path to the outgoing residual. The signed TV witness also satisfies the required assumptions.

Specifically choose \(s_c=\operatorname{sign}(a_c-b_c)\) from the projected results and hold it fixed. Set

\[
F_j=\tfrac12\big[s_{c(j)}-(K^Ts)_j\big],
\qquad \|F\|\le1.
\]

Its projected expectation is exactly the projected TV. Outside the declared retained support the kernel is the identity and \(F=0\). For the exact full-law state,

\[
D_{\rm TV}^{\rm exact}\ge\langle F\rangle_{\rm exact}
\ge D_{\rm TV}^{P}-(2ed+e^2).
\]

The resulting lower bound is valid for this explicitly specified reference even though the witness was chosen after inspecting the calculation. Its expectation's upper endpoint is not a refined upper bound on TV; that distinction is retained in the output.

All truncation, contraction and dual inequalities are analytic. Numerical exponentials, norms and eigendecompositions use ordinary floating point. This packet does not provide an interval-arithmetic certificate for floating-point roundoff. Published decimal intervals are rounded outward.

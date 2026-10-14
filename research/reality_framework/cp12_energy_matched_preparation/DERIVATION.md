# Derivation: phase-only energy matching and full-law initial response

All statements use the frozen CP-12 couplings unless a symbol is retained. Put M=binom(N,2), d=Q(M−Q), D=binom(M,Q). Let P be the neutral sea subspace, one state for each Q-edge graph. On P, the diagonal mass/electric/bond energy is κQ. Hopping maps P into its orthogonal charged sector. The rewiring term preserves P and equals the Johnson-graph Laplacian

\[
R=\frac{\Gamma}{M}(dI-A_J).
\]

A_J connects Q-subsets differing by one removed and one inserted connection. The full matter Hamiltonian does not possess all permutations of the M potential connections. This larger symmetry is used only inside R and for the specifically justified initial moments.

## 1. Mean energy and variance

Let z_e=e^{iθ} for e in S and z_e=1 otherwise, and ψ(G)=D^−1/2 product_{e in G} z_e. For e≠f the probability, in the uniform Q-subset measure, that e is occupied and f absent is Q(M−Q)/[M(M−1)]. Therefore

\[
\langle A_J\rangle=\frac{d}{M(M-1)}\sum_{e\ne f}z_f\overline{z_e}
=\frac{d}{M(M-1)}\left(\left|\sum_ez_e\right|^2-M\right).
\]

Since |sum z_e|²=M²−2a(M−a)(1−cosθ), imposing |sum z_e|²=M gives the matching angle in REPORT.md and <R>=Γd/M. This is also the incoherent mixture's rewiring energy. Hopping has zero expectation in every neutral preparation.

For any normalized neutral vector ψ,

\[
\operatorname{Var}_\psi H=\operatorname{Var}_\psi R+2t_{\rm hop}^2Q.
\]

The charged states created by the 2Q individual hops on any fixed graph are orthogonal, including their electric fields; different graphs remain orthogonal. Hence ||T_hop ψ||²=2t_hop²Q independent of phases. Neutral and charged contributions to (H−<H>)ψ are orthogonal. Any two masks of the same cardinality are related by a permutation of the M connection coordinates, which preserves R. Their R moments, and therefore full-H mean and variance, agree. No claim about all moments of full H follows.

For the dephased mixture,

\[
\operatorname{Var}_{\rm mix}H=2t_{\rm hop}^2Q+\Gamma^2d/M^2.
\]

## 2. Why full-H graph derivatives reduce through order two

Let O be graph-diagonal, acting as a scalar on each graph's matter/electric fibre. It commutes with hopping and the non-rewiring diagonal terms. Thus [H,O]=[R,O]. For neutral initial states, [R,O] preserves P; one-hop charged states are orthogonal to P. The cross terms involving T_hop in <[H,[H,O]]> vanish. The remaining diagonal part acts as κQ on P and also drops out. Consequently

\[
\langle O(T)\rangle=\langle O\rangle
+iT\langle[R,O]\rangle-\frac{T^2}{2}\langle[R,[R,O]]\rangle+O_N(T^3).
\]

This holds at frozen hopping t_hop=1 with unbounded electric fields. It is not a weak-hopping or Born-Oppenheimer approximation. Neutral sea vectors have finite support and finite moments at each fixed N. The remainder is fixed-N; no thermodynamic interchange of limits is asserted. At higher order, matter/electric terms generally survive, as the exact T⁴ certificates demonstrate.

## 3. Two-color quotient

Write b=M−a, k=|G intersect S|. The number of such graphs is C(a,k)C(b,Q−k). The normalized uniform state in each class gives an orthonormal basis |k>. Put

\[
p_k=(Q-k)(a-k),\qquad m_k=k(b-Q+k).
\]

They count transitions from k to k+1 and k−1. Removing transitions within the same class leaves

\[
R_{kk}=\frac{\Gamma}{M}(p_k+m_k),\quad
R_{k,k+1}=-\frac{\Gamma}{M}\sqrt{(k+1)(a-k)(Q-k)(b-Q+k+1)}.
\]

The initial quotient vector has components sqrt(w_k)e^{iθk}, with
w_k=C(a,k)C(b,Q−k)/C(M,Q). A global phase is immaterial. Numerics use mode-centered recurrence for w_k and a recurrence for phases; subtracting large log-gamma values or exponentiating large θk separately introduces avoidable cancellation in tiny response coefficients.

If v₁=Rψ and v₂=R²ψ, then the k-probability coefficients are

\[
p_k^{(1)}=2\operatorname{Im}(\overline\psi_kv_{1,k}),\qquad
p_k^{(2)}=|v_{1,k}|^2-\operatorname{Re}(\overline\psi_kv_{2,k}).
\]

A constant energy shift may be removed in this evaluation. The quotient is exact for R at all times but is used for full H only as justified above.

## 4. Independent rational formula and the degree-correlation limitation

The eigenvector amplitudes need not be numerically represented. Put g=Γ/M, c=cosθ and s²=1−c². For a single graph of color count k,

\[
(R\psi)(G)/\psi(G)=A_k+i\sin\theta B_k,
\quad A_k=g(p_k+m_k)(1-c),\quad B_k=-g(p_k-m_k).
\]

Thus its second-order relative probability coefficient is the polynomial

\[
q_2(k)=A_k^2+s^2B_k^2-g\{(p_k+m_k)A_k
-p_k(cA_{k+1}-s^2B_{k+1})
-m_k(cA_{k-1}+s^2B_{k-1})\}.
\]

Although written as a quartic expression, cancellations make it quadratic. For balanced pools a=b=M/2 and the matching angle c=−1+2/M,

\[
\boxed{q_2(k)=\frac{4\Gamma^2(M-1)^2}{M^4}
\left[4(M-1)(k-Q/2)^2-Q(M-Q)\right].}
\]

Its uniform mean is zero since Var(k)=Q(M−Q)/[4(M−1)]. All main Q=3N sizes have balanced pools.

For a graph observable invariant under event permutations, E[O n_e] is independent of the edge e. For two distinct potential edges, E[O n_en_f] has only two types: incident and disjoint. Since k²=k+2 sum_{e<f in S} n_en_f, E[O q₂(k)] depends on the mask only through its number of incident marked-edge pairs,

\[
W(S)=\sum_i\binom{d_i(S)}2.
\]

The first-order coefficient vanishes by the one-edge symmetry. Therefore the initial invariant graph response at order T² cannot encode more than this mask degree statistic. Different short-cycle coefficients at this order do not constitute independent macroscopic topology control. This argument is for fixed-size masks and the stated preparation family, not arbitrary coherent graph states.

For an r-edge motif with n_{r,j}(S) embeddings using j marked connections, the conditional expectation is

\[
O_r(k)=\sum_{j=0}^r n_{r,j}(S)
\frac{(k)_j}{(a)_j}\frac{(Q-k)_{r-j}}{(b)_{r-j}}.
\]

Here (x)_j is a falling factorial. Both O_r and q₂O_r are polynomials. Their averages are evaluated exactly with

\[
E[(k)_j]=(Q)_j(a)_j/(M)_j,
\quad P(k)=\sum_j\frac{\Delta^jP(0)}{j!}(k)_j.
\]

This rational route is independent of the floating-point tridiagonal quotient. The checker also verifies that the finite-difference sequence terminates at the specified degree, counts colored four-cycles by explicit enumeration for N=4,6,8, and compares the balanced closed expression with the quotient coefficient. Results store exact rational coefficients as well as decimals.

## 5. Exact small-system full-H certificate

For N=4,Q=4, z=−2/3+i sqrt(5)/3. Every initial amplitude, apart from the common 1/sqrt(15), belongs to Q(i sqrt(5)). A pair (a,b) represents a+i sqrt(5)b; multiplication is
(a,b)(c,d)=(ac−5bd,ad+bc).

The independent published `exact_destination.transitions` routine gives integer coefficients for a common scale lcm(8,10M)=120. Apply it four times to all 15 neutral sea amplitudes, retaining every reachable physical state. For even j,

\[
c_j=\frac1{15}\sum_{a=0}^j\frac{i^a(-i)^{j-a}}{a!(j-a)!}
\langle H^a\psi_{\rm unnorm}|O|H^{j-a}\psi_{\rm unnorm}\rangle.
\]

The common denominator scale^j is restored in code. Odd coefficients are represented as a rational multiple of sqrt(5) and vanish for j=1,3 in the examples. This independently verifies the full-law, unregulated energy moments and square coefficients, without importing the regulated matrix builder or using the two-color quotient.

## 6. Information and physical scope

For every pure preparation here, H∞=0, S(ρ)=0, and the diagonal graph distribution has entropy log₂D. Hence neutral-basis relative-entropy coherence S(dephase ρ)−S(ρ)=log₂D for all masks. The two states are indistinguishable by these scalar resources but are physically distinct relative to the unchanged H.

The gate is a preparation intervention, not a term secretly added to CP-12 evolution. It is not an endogenous state-selection principle. There is no preparation from the flat mixture by this unitary, and no derivation of a cosmological initial condition. The high-density study computes all-square counts and degree wedges, not extensive vertex-disjoint loop packing, intrinsic three-dimensional volume growth, helicities, or propagation cones.

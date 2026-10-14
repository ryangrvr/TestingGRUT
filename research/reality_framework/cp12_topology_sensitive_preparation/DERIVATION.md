# Exact cubic transfer and topology discrimination

## 1. Full-law third-order reduction

Let P project onto neutral seas, one vector for every Q-edge graph. Subtract the constant κQ from H for this calculation; observable evolution is unchanged. Write H=A+T+R_full, where A is the positive mass/electric diagonal, T the hopping operator with strength t_hop, and R_full the positive rewiring Laplacian. On P, A=0 and R_full=R, the neutral Johnson Laplacian.

For an initial ψ in P, Tψ consists of one vertex fermion, one edge hole and one nonzero electric half-leg. Let K₁ be this one-hop sector. Its A energy is Δ=2m+ε_E/2 and its rewiring block is R₁. Define

\[
\alpha=2Q t_{\rm hop}^2,\qquad \beta=2(Q-1)t_{\rm hop}^2.
\]

The exact identities are

\[
PT^2P=\alpha P,\qquad PTR_1TP=\beta R.
\]

**Counting proof.** Each neutral graph has 2Q orthogonal single-hop outputs. A second hop returns to a neutral sea only by reversing the same hop, including its fermionic sign; thus PT²P=2Q t²P. In each single-hop state the excited edge is pinned: its center is empty and one electric leg is nonzero. The other Q−1 edges are removable. For a neutral graph swap G→G−e+f, there are 2(Q−1) common hopping channels on its unchanged edges. The fermionic signs before and after that swap agree because the sea occupations and the pinned hop coincide. This gives β times the neutral off-diagonal rewiring entry. For the diagonal, summing the 2Q single-hop channels, each with (Q−1)(M−Q) swaps, gives exactly β times R's neutral diagonal. No destination restriction or graph diagonal is omitted.

Put v_j=H^jψ. Their neutral components are

\[
Pv_1=R\psi,\quad Pv_2=(R^2+\alpha)\psi,
\quad Pv_3=[R^3+(2\alpha+\beta)R+\alpha\Delta]\psi.
\]

The one-hop component of v₂ is (Δ+R₁)Tψ+TRψ. Its other component has two vertex excitations and cannot overlap with v₁. This remains true with unbounded rotors: the argument is by occupation sectors, not an electric cutoff.

Let O be real graph-diagonal, constant on each graph's matter/electric fibre. It commutes with T and A. Define r=<ψ|OR|ψ>. The Taylor coefficients are

\[
c_1^H=2\operatorname{Im}r,
\quad c_3^H=\operatorname{Im}\langle v_1|O|v_2\rangle
-\frac13\operatorname{Im}\langle\psi|O|v_3\rangle.
\]

The hopping terms in the first imaginary part sum to β Im r: the two α terms are conjugates and cancel, and αΔ is real. In the second imaginary part they sum to (2α+β) Im r. Hence

\[
c_3^H=c_3^R+\frac23(\beta-\alpha)\operatorname{Im}r
=c_3^R-\frac23t_{\rm hop}^2c_1^R.
\]

The parent's flat-amplitude preparations have ψ(G)=D^−1/2 z^{k_S(G)}. Their first graph probability correction is linear in k. For permutation-invariant O, E[O n_e] is independent of the potential edge, and sum_e n_e=Q is fixed. Therefore c₁(O)=0. Together with the parent's order-two proof, this establishes full-neutral agreement through T³ for all such O in this preparation family. It makes no finite-time substitution claim.

The independent full-moment checker also tests non-invariant O=k_S. Its nonzero correction is −4√15/375, exactly as the identity requires.

## 2. Exact cubic probability polynomial

Let a=|S|, b=M−a, g=Γ/M, c=cosθ and s=sinθ. Define p_k=(Q−k)(a−k), m_k=k(b−Q+k). For an unnormalized graph amplitude ψ(G)=z^k, write (R^jψ)(G)/ψ(G)=u_j(k), u₀=1. Then

\[
u_{j+1}(k)=g[(p_k+m_k)u_j(k)-p_k z u_j(k+1)-m_k\bar z u_j(k-1)].
\]

Represent u_j=A_j+i s B_j, so all coefficients are rational functions of M,Q,a,c. The cubic probability correction is

\[
q_3(k)=\operatorname{Im}[\overline{u_1(k)}u_2(k)]
-\frac13\operatorname{Im}u_3(k).
\]

Formal expansion cancels all powers above k³. `symbolic_check.py` establishes the identity for formal variables, not only selected numerical parameters. The k³ coefficient in q₃/s is

\[
b_3=\frac43g^3(1-c)
[8(1-3c)-2M(1-13c)-3M^2(1+3c)+M^3(1+c)].
\]

It depends on Q and a only through the supplied matching c. This is a symbolic polynomial identity, not a fitted effective interaction.

## 3. Which three-edge statistics differ

The main masks use copies of the parent's six-event pair plus identical complete bipartite cross-block connections. They share the degree of every vertex and the statistic

\[
J(S)=\sum_{\{u,v\}\in S}(d_u-1)(d_v-1).
\]

For internal edges, both six-event masks have the same number of endpoint-degree types. Every vertex in a block receives the same cross-block degree shift. All cross-block edges also have the same endpoint degrees in both masks. Therefore J agrees exactly, not just numerically.

Let Δτ be triangle-rich-minus-square-rich marked triangle count. Among unordered subsets of three marked edges, the changes in counts of triangle, three-star, three-edge path, wedge-plus-disjoint-edge, and three-disjoint-edge patterns are respectively

\[
(\Delta\tau,\ 0,\ -3\Delta\tau,\ 3\Delta\tau,\ -\Delta\tau).
\]

Proof: the three-star count is sum_v binom(d_v,3), fixed. J equals the number of simple three-edge paths plus 3τ, so the path change is −3Δτ. W(S)(a−2) counts a marked incident pair and a third edge; a triangle or three-star contributes three, a path two and a wedge-plus-edge one. This fixes its change at +3Δτ. Fixed total binom(a,3) then fixes the matching change at −Δτ.

The block family has Δτ=2(N/6)=N/3. Its marked-square difference is −N/3, but marked four-edge topology is not independently accessible to q₃.

## 4. Closed overlap formulas

Put p_r=(Q)_r/(M)_r, the probability that r fixed potential connections are all occupied in a uniform Q-edge graph. Lower-order k moments cancel between the masks. Only the ordered three-marked-edge term in k³ can distinguish them, giving a factor 3!=6.

For the combination of selected three-edge patterns

\[
\text{triangle}-3\text{path}+3(\text{wedge+edge})-\text{matching},
\]

count overlap with an observed triangle. The combined numbers of observed triangles sharing 3,2,1,0 selected edges are (1,−3,3,−1). Thus

\[
\Delta E[C_3 k^3]=6\Delta\tau[p_3-3p_4+3p_5-p_6]
=6\Delta\tau\frac{(Q)_3(M-Q)_3}{(M)_6}.
\]

The last equality is inclusion-exclusion for three specified occupied and three specified absent edges. It requires no independence approximation.

For observed squares, a selected three-edge path closes in one way, whereas the other selected patterns cannot supply three edges of a square. The combined overlap count for 3 shared edges is −3. A fixed incident selected edge pair belongs to N−3 squares, and a disjoint pair to two. Since the incident-pair count cancels, the 2-shared-edge count is +9. Fixed single-edge incidence then gives −9 for one shared edge; fixed total square count gives +3 for zero. Consequently

\[
\Delta E[C_4 k^3]=-18\Delta\tau[p_4-3p_5+3p_6-p_7]
=-18\Delta\tau\frac{(Q)_4(M-Q)_3}{(M)_7}.
\]

Multiplying by s b₃ gives the two closed response formulas in REPORT.md. Their ratio is −3(Q−3)/(M−6). The two-edge degree-wedge observable has zero change under the same three-pattern combination, as confirmed independently in all calculations.

An equivalent calculation conditions on k using colored motif counts n_{r,j}(S):

\[
E[C_r|k]=\sum_j n_{r,j}(S)\frac{(k)_j}{(a)_j}\frac{(Q-k)_{r-j}}{(b)_{r-j}}.
\]

Average products with q₃ using E[(k)_j]=(Q)_j(a)_j/(M)_j. The code implements this separate route with exact rational factorial moments and checks equality with the overlap formulas.

## 5. Controlled asymptotic statement

For the declared family Q=3N, M=N(N−1)/2 and a=M/2+O(1), energy matching gives

\[
c=-1+2/M+O(M^{-2}),\quad s\sim2/\sqrt M,
\quad b_3\sim\frac{64}{3}\frac{\Gamma^3}{M}.
\]

Also 6(Q)_3(M−Q)_3/(M)_6∼1296 N^−3 and Δτ=N/3. Hence

\[
\Delta c_3(C_3)\sim36864\sqrt2\Gamma^3 N^{-5}.
\]

The exact square ratio tends to −18/N, giving the N^−6 formula. This is a derivative asymptotic along the declared family. It does not bound higher derivatives uniformly and must not be extrapolated to times proportional to N^{5/3}, or any other growing timescale, without remainder control.

## 6. Phase reversal and physical boundary

For real H and O in the frozen basis, conjugating the initial vector sends <O(T)> to <O(−T)>. Writing conjugation as the antiunitary K, K H K=H implies K 1_B(H) K=1_B(H) for every real spectral projector. Thus the conjugate preparations have identical energy spectral measures, not only identical low-order moments. Initial graph probabilities and scalar information resources are unchanged; every odd Taylor coefficient reverses sign. This full-energy-distribution statement is for each state and its conjugate; it is not asserted for the two distinct-topology masks. No irreversible assembly, state-selection principle or time arrow follows from the cubic coefficient.

The controller mask and coherent neutral input are external preparations. The gate exp(iθ sum_{e in S}n_e) preserves Q and Gauss's law, but does not purify a flat mixture or select its own S. Intrinsic dimension, disjoint loop packing and local gauge-invariant source-detector propagation have not been calculated here. All theorem transfers require the stated neutral initial vectors and graph-only observables.

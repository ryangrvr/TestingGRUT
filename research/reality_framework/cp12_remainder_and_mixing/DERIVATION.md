# CP-12: full-Hamiltonian bounds on regular geometry and its mixing

This is an additive analysis of frozen CP-12, based on
`e7a6535e849d876a0f7b4611c0a2d98ae6437148`. No interaction, coupling,
state-counting rule, or canonical research card is changed. The bounds below
refer to the Hamiltonian in `research/constructive_physics/attempt12/model.py`:

\[
m=t=1,\quad \kappa=\tfrac12,\quad \epsilon_E=\tfrac14,
\quad\Gamma=\tfrac1{10},\quad M=\binom N2.
\]

Electric fields are integers on the two half-legs of each potential edge.
Absent edges have filled edge fermions and zero fields. Present edges obey
the edge Gauss law; vertices obey their own Gauss law. There are exactly
\(Q\) present edges and \(M\) fermions. The energy convention includes the
frozen \(\kappa Q\) term; it is not the neutral-subtracted cluster energy.

At finite \(N\), the electric quadratic form is confining, and the hopping
and graph operators are bounded. The unregulated physical Hamiltonian has
a ground state. Finite electric cutoffs are used only for the executable
matrix controls, not for the analytical phase bounds.

## 1. Exact rewiring fibers

Fix the fermion occupations and electric fields. An edge is *locked* if its
edge fermion is empty or either electric half-leg is nonzero. Write \(S\)
for the number of such edges. They must all be present. The remaining
\(K=M-S\) eligible edge positions contain \(r=Q-S\) neutral present bonds.
Their rewiring graph is exactly the Johnson graph \(J(K,r)\): one removes
one occupied position and adds one unoccupied position.

The positive rewiring operator on this fiber is

\[
R_f=\frac\Gamma M\mathcal L[J(K,r)].
\]

Its eigenvalues are \(\Gamma j(K+1-j)/M\), for
\(0\le j\le\min(r,K-r)\). The zero eigenvector is the uniform superposition
of compatible graphs. Its first nonzero eigenvalue is \(\Gamma K/M\).
Let \(P\) project onto the direct sum of these uniform fiber vectors. For
\(0<Q<M\),

\[
R\ge\frac{\Gamma(M-Q+1)}M(I-P).
\]

This is a gap of the **rewiring operator**, not a gap of the full Hamiltonian.
The Johnson/token-graph spectrum is established mathematics; see Reyes,
Dalfó and Fiol, [arXiv:2310.16929](https://arxiv.org/abs/2310.16929), in
particular the complete-graph token Laplacian's first gap. The application
here follows from CP-12's unchanged neutral-removal rule.

The entire graph-changing off-diagonal operator \(V\) has the exact norm

\[
\boxed{\|V\|=\Gamma Q(M-Q)/M.}
\]

Each fiber has at most \((Q-S)(M-Q)\) neighbors. The completely neutral
fiber attains the maximum with its uniform vector. Thus the incoming
maximum-norm bound is sharp as an operator statement, although it is not
an estimate of an actual ground-state energy shift.

There is also an exact projected Hamiltonian \(PHP\). Its diagonal is
\(\kappa Q+2mN_v+(\epsilon_E/2)\sum E^2\). Creating a new locked edge
multiplies the frozen hopping matrix element by
\(\sqrt{(Q-S)/(M-S)}\); removing a lock multiplies it by
\(\sqrt{(Q-S+1)/(M-S+1)}\). Hops leaving the lock count unchanged have
factor one. These factors count compatible neutral graphs. Projection is
an exact compression, not an assertion that the full ground state is in
this kernel.

## 2. A sharper regular-sector coupling bound

Let \(P_d\) include every exactly \(d\)-regular graph and every allowed
matter/electric state on it, at \(Q=dN/2\). One edge replacement always
changes some degrees, so \(P_dVP_d=0\). Its compressed Hamiltonian is a
direct sum of complete graph blocks, **including** their rewiring diagonal.

Any regular state has at most \(Q(M-Q)\) outgoing replacements. A
nonregular graph reached in one step has at most \(d+1\) regular
predecessors:

* Disjoint removed/added edges leave two deficit and two excess vertices.
  Their reverse edges are uniquely specified.
* Edges sharing one endpoint leave one deficit and one excess vertex.
  A reverse operation chooses a neighbor of the excess vertex, which has
  degree \(d+1\).

Locking and Gauss constraints can only delete transitions. The row-column
Schur bound therefore gives

\[
\boxed{b_N:=\|P_dH(I-P_d)\|
 \le\frac\Gamma M\sqrt{(d+1)Q(M-Q)}=O(N^{-1/2}).}
\]

In fact, its norm equals the corresponding neutral graph matrix's norm:
every physical fiber gives a submatrix of that matrix, and the neutral
fiber itself is included. This controls the full physical space, including
arbitrarily large integer fields.

If \(E_R\) and \(E_C\) are the bottoms of the regular and complementary
compressed sectors, then

\[
\min(E_R,E_C)-b_N\le E_0\le\min(E_R,E_C).
\]

Internal mixing in the complementary sector must still be included. It
cannot be replaced by its best individual graph energy.

For any initial state supported in either sector, including a coherent
superposition, transition probability to the other is at most \(b_N^2T^2\).
If norm-one source and detector operators commute with \(P_d\), and the
initial state is in the regular sector, then

\[
|R_{\rm full}(T)-R_{\rm regular\ block}(T)|\le4b_N^2T^2.
\]

For the last statement, write the compressed propagator as \(P_dUP_d\).
Its difference from evolution by \(P_dHP_d\) is at most \(b_N^2T^2/2\);
the complementary amplitude is at most \(b_NT\). The compressed detector
error is at most \(2b_N^2T^2\), and taking its commutator with the source
adds another factor two. This protects a **prepared** coherent regular
ensemble. Conversely, assembly into exact regularity from a state wholly
outside it requires \(T=\Omega(\sqrt N)\) for fixed nonzero success
probability. Neither conclusion defines a distance or a physical clock.

## 3. A Gauss-law lower bound for every graph block

The true diagonal coefficient per removable edge is

\[
a_N=\Gamma\frac{M-Q}{M}.
\]

It is not an arbitrary offset: the missing off-diagonal transitions leave
their positive degree diagonal in each graph compression.

Split a graph block into vertex stars. Star \(i\) receives half the cost of
each incident edge, all its vertex mass, half each incident edge's electric
energy, half its removable-edge diagonal, and all hopping to that vertex.
Summing the stars exactly reconstructs the graph block.

Edge Gauss law reads \(E_{ie}+E_{je}=1-n_{b_e}\). For
\(a_N\le\epsilon_E/2\), every edge's contribution to a star is at least
\(a_N/2\). If \(n_{a_i}=1\), vertex Gauss law requires
\(\sum_{e\ni i}E_{ie}=1\). At least one own half-leg is nonzero, adding
at least \(\delta=\epsilon_E/4-a_N/2\). Hence the star's diagonal is at
least

\[
\frac{(\kappa+a_N)d_i}{2}+(2m+\delta)n_{a_i}.
\]

For the combined incident fermion
\(F_i=\sum_{e\ni i}\mathrm{orient}_{ie}U_{ie}b_e\), the canonical
anticommutator is \(\{F_i,F_i^\dagger\}=d_i\). The off-diagonal coupling
between the two central occupation sectors has norm at most \(t\sqrt{d_i}\).
The elementary two-block inequality then gives

\[
h_i\ge f_a(d_i),\qquad
f_a(d)=\frac{(\kappa+a)d}{2}
 +\frac{\Delta_s-\sqrt{\Delta_s^2+4dt^2}}{2},
\quad\Delta_s=2m+\epsilon_E/4-a/2.
\]

The function \(f_a\) is convex for \(d\ge0\). Thus at \(Q=2N\), Jensen's
inequality gives **every** complete fixed-graph block, regular or otherwise,
the lower bound

\[
H_G\ge N\ell_N,\qquad \ell_N=f_{a_N}(4).
\]

For the frozen constants and \(N\to\infty\),

\[
\boxed{\ell_\infty=\frac{353-\sqrt{128321}}{160}
 =-0.03262004144948\ldots .}
\]

This covers all hopping paths and all electric fields. It uses no
weak-hopping expansion. `star_bound_check.py` reconstructs the individual
stars from complete regulated physical matrices and checks the diagonal
inequality, coupling norm, and exact sum. Its degree-four central-star
matrix attains the local lower bound to floating-point accuracy.

## 4. A coherent variational state below that lower bound

Start from the neutral fermion sea. A *matching* selects \(k\) distinct
vertices, each paired with a distinct incident edge. Move that edge fermion
to its vertex and put electric field \(+1\) on that half-leg. All other
fields vanish. Every configuration satisfies both Gauss laws exactly.

For each matching, take the uniform superposition of all \(Q\)-edge graphs
containing its \(k\) locked edges. This lies exactly in \(\ker R\). Different
matchings have orthogonal matter/field configurations. The phases are
fixed by the frozen hopping orientation and fermionic signs; disjoint
even fermion bilinears commute. No tunable phase or new interaction is
introduced.

The number of matchings is

\[
m_k=\binom Nk H_k,\quad H_0=1,\ H_1=N-1,\quad
H_k=(N-1)H_{k-1}-(k-1)H_{k-2}.
\]

To derive it, initially let every selected vertex choose one of \(N-1\)
incident edges. The only collision is two endpoints choosing their common
edge. Inclusion-exclusion over disjoint pairs gives

\[
H_k=\sum_{j=0}^{\lfloor k/2\rfloor}
(-1)^j\frac{k!}{2^j j!(k-2j)!}(N-1)^{k-2j}.
\]

The normalized uniform matching states \(|u_k\rangle\) have a tridiagonal
compression of the **full** Hamiltonian:

\[
T_{kk}=\kappa Q+\Delta k,\qquad \Delta=2m+\epsilon_E/2=17/8,
\]
\[
\boxed{T_{k,k+1}=t(k+1)
\sqrt{\frac{Q-k}{M-k}\frac{m_{k+1}}{m_k}}.}
\]

Every \((k+1)\)-matching has \(k+1\) parents; the square root also contains
the compatible-graph occupation factor from section 1. Other hopping
channels leave this variational subspace. It need not be invariant:
Rayleigh-Ritz still gives an upper bound for the complete ground energy.
The projection and \(R|u_k\rangle=0\) are independently checked against
full physical matrices with 30, 872, and 3,975 states.

For \(Q=\nu N\), \(k=xN+o(N)\),

\[
T_{k,k+1}/N\longrightarrow t\sqrt{2x(1-x)(\nu-x)}.
\]

Indeed \(H_{k+1}/H_k=(N-1)-k/(H_k/H_{k-1})\), and for \(N\ge5\),
\(N-3\le H_k/H_{k-1}\le N-1\), by induction. This also supplies a
conservative finite-size coupling bound without enormous integer counts.

Choose an alternating-sign tent wavepacket centered at \(k=N/4\) with
width \(\lfloor\sqrt N\rfloor\). Adjacent coefficient overlap tends to one;
its support concentrates at \(x=1/4\). At \(\nu=2\), its Rayleigh quotient
\(U_N\) obeys

\[
\boxed{\lim_{N\to\infty}U_N/N=u_\infty
=49/32-\sqrt{21/8}=-0.08893517460197\ldots .}
\]

Thus \(E_0\le U_N\) and

\[
g_\infty:=\ell_\infty-u_\infty=0.05631513315248\ldots>0.
\]

This trial is an allowed state, not a proposed preparation or identification
of the actual ground state. Its extensive energy advantage establishes
that graph coherence matters at the frozen coupling, beyond the tiny
fixed-graph cluster comparisons.

## 5. Exact regular-ground weight tends to zero

Section 2 and the extensive gap just obtained imply for any normalized
ground state \(|\psi_0\rangle\),

\[
(P_4HP_4-E_0)P_4\psi_0=-P_4H(I-P_4)\psi_0.
\]

If \(G_N=N\ell_N-U_N>0\),

\[
p_4:=\|P_4\psi_0\|^2\le
\frac{b_N^2}{G_N^2+b_N^2}\le\frac{b_N^2}{G_N^2}.
\]

Here \(b_N^2\sim0.2/N\), while \(G_N\sim g_\infty N\). Consequently

\[
\boxed{p_4=O(N^{-3})\longrightarrow0.}
\]

The finite certificates use exact integer matching counts, rational squared
couplings, and 80-digit outward Decimal bounds on square roots, Rayleigh
quotients, and gaps. The optional floating-point tridiagonal eigenvalue is
only a diagnostic, not the quantity used to certify the bound.

## 6. Extension to nearly regular graphs

Exact regularity alone is too restrictive a physical target. Define the
absolute degree defect

\[
D(G)=\sum_i|d_i(G)-4|,
\]

and let \(P_K\) include all graphs with \(D\le K\), with their complete
physical matter/electric spaces. A swap changes \(D\) by at most four.
There are at most \(K\) defective vertices, whose summed degrees are at
most \(5K\).

At any graph in this sector, swaps touching an already defective vertex
number no more than

\[
C_N(K)=5K(M-Q)+QK(N-1).
\]

This counts a defective endpoint in the removed edge or in the added edge;
overcounting only strengthens the upper bound. All other swaps increase
\(D\) by two or four.

Give each graph the positive Schur weight \(\eta^{D(G)}\), \(0<\eta\le1\).
The nonnegative adjacency's weighted row sums are at most

\[
\eta^2 Q(M-Q)+\eta^{-4}C_N(K).
\]

This proves the internal off-diagonal norm bound

\[
\|P_KVP_K\|\le r_N(K):=
\min\left\{\frac{\Gamma Q(M-Q)}M,
\inf_{0<\eta\le1}\frac\Gamma M
[\eta^2Q(M-Q)+\eta^{-4}C_N(K)]\right\}.
\]

For \(K=o(N)\), choosing \(\eta\) near \((9K/N)^{1/6}\) gives
\(r_N(K)=O(N^{2/3}K^{1/3})=o(N)\). The internal sector's energy is thus
at least \(N\ell_N-r_N(K)\).

A graph outside this sector with a predecessor inside has \(D\le K+4\).
Any reverse swap lowering \(D\) touches one of its defective vertices.
Applying the ordinary row-column bound across the sector boundary yields

\[
\boxed{\|P_KH(I-P_K)\|^2\le
\frac{\Gamma^2}{M^2}Q(M-Q)C_N(K+4)
=O(N(K+4)).}
\]

These graph bounds apply to the physical fibers as submatrix bounds.
Since \(N\ell_N-r_N(K)-U_N\sim g_\infty N\), the projected ground-state
equation gives

\[
\boxed{\Pr_{\psi_0}[D\le K_N]
=O((K_N+4)/N)\longrightarrow0\quad\text{for }K_N=o(N).}
\]

Thus a subextensive number of absolute degree defects cannot dominate the
ground-state graph distribution. This is stronger than the exact-regular
result, though it remains a statement about approach to degree four.

One quantitative asymptotic corollary is possible. At \(K=\lfloor\delta N\rfloor\),

\[
\limsup r_N(K)/N\le\beta\delta^{1/3},\quad
\beta=3\Gamma\,9^{1/3}=0.6240251469\ldots,
\]
\[
\limsup \|P_KH(I-P_K)\|^2/N^2\le18\Gamma^2\delta=0.18\delta.
\]

Whenever \(g_\delta=g_\infty-\beta\delta^{1/3}>0\),

\[
\limsup\Pr[D\le\lfloor\delta N\rfloor]\le
\frac{0.18\delta}{g_\delta^2+0.18\delta}.
\]

For \(\delta=10^{-4}\), conservative rational comparisons give
\(g_\infty>0.056\), \(\beta<0.625\),
\(\delta^{1/3}<0.047\), and hence the last probability is below \(1/40\).
It follows that

\[
\boxed{\liminf_{N\to\infty}\langle D\rangle/N
\ge39/400000=9.75\times10^{-5}.}
\]

This deliberately weak positive density bound is rigorous within the
frozen model; it is not an experimentally inferred defect fraction.
`near_regular_exclusion.py` checks the weighted bounds against all 455
graphs at \(N=6,Q=12\), verifies the exact rational comparisons, and
produces conservative finite certificates up to \(N=2^{24}\).

## 7. What the finite remainder controls actually establish

The complete periodic graph blocks \(C_4,C_6,C_8\), at the inherited host
coefficient \(a=0.1\), are independently diagonalized at electric cutoffs
3 and 4. Their complete finite-host linked sums reproduce those energies.
The edge-seven remainder on \(C_8\) is approximately \(8.16\times10^{-8}\)
per vertex in the neutral-subtracted convention. It is a measured finite
host remainder, not an upper bound on the infinite square-lattice tail.

For graph mixing, split the full finite Hamiltonian into \(P=\ker R\) and
its complement:

\[
H=\begin{pmatrix}A&B^\dagger\\B&C\end{pmatrix},\qquad
F(z)=A-B^\dagger(C-z)^{-1}B-zI.
\]

For \(z<\inf\operatorname{spec}C\), a ground eigenvalue is a zero of the
lowest eigenvalue of \(F\). If \(\operatorname{spec}(C-z)\subset[a,b]\),
\(a>0\), the degree-\(k\) Chebyshev approximation to its inverse has error

\[
\epsilon_k\le\frac{2}{\sqrt{ab}}\frac{r^{k+1}}{1-r},
\quad r=\frac{\sqrt b-\sqrt a}{\sqrt b+\sqrt a}.
\]

This follows by summing the geometric tail of the Chebyshev series
\(1/\sqrt{ab}[1+2\sum_{j\ge1}(-r)^jT_j]\). The Schur norm error is
at most \(\|B\|^2\epsilon_k\). Moreover
\(F'(z)=-I-B^\dagger(C-z)^{-2}B\le-I\), so approximation error plus
root defect transfers to an energy enclosure, provided it lies below
the complementary spectrum.

The complete \(N=4,Q=3,L=3\) test uses all 1,032 physical states and
recovers the direct ground energy within about \(10^{-14}\). At degree
128 its analytic polynomial error is below \(4\times10^{-14}\); the
reported energy enclosure retains a \(10^{-10}\) floating-point floor.
Its spectral enclosure uses residual-padded dense eigenvalues, **not**
interval arithmetic. This is a validation of a conditional remainder
method, not a large-system certified spectral enclosure. None of these
finite numerical controls is needed to truncate the phase proof above.

## 8. Scientific scope

The zero-temperature phase result concerns frozen CP-12 at fixed
\(Q/N=2\), its supplied all-pairs \(\Gamma/M\) rewiring, and its specified
physical Hilbert space. It excludes ground-state concentration on exact
four-regular graphs and on graphs whose total absolute degree defect is
subextensive. It includes balanced lattices, thin tori, regular expanders,
and regular fragments without needing to rank them.

It does **not** exclude a spatial phase with a nonzero density of degree
defects, a Poisson-like manifold discretization, a metastable prepared
geometry, other bond-density sectors, or every emergent-spacetime model.
The trial state does not identify the actual alternative phase, and the
analysis derives neither a Lorentzian metric nor Einstein gravity. The
same dynamics can protect a prepared regular ensemble while failing to
select it as its ground phase. Those are different physical requirements.

# Extensive circulation: all-time probability bounds and charged controls

Additive continuation of `36138c6140062b20c60864bce12a893c91dca8eb`. Frozen CP-12, Q conservation, Gauss constraints and all couplings are unchanged. No theorem below establishes formation or a causal cone. Numerical checks validate analytic steps; the all-size statements rest on the proofs.

## 1. Setup and target sectors

Write H=kappa Q+D+R+K, where D=2m n_v+F, F=(epsilon_E/2)sum E^2, R is the positive neutral-rewiring Laplacian, and K=t(A+A^dagger) is hopping. A creates one vertex fermion and one edge-mode hole. At fixed finite N, K and R are bounded; only the electric diagonal is unbounded.

Let mathcal C be the union of directed electric-circulation cycles, e_C its edge count, and h_C its edge-hole count. The parent proves the diagonal inequality

\[
\sum E^2\ge2e_C+n_v+2h_C.
\]

Choose a configuration event G determined by circulation/fermion variables, implying e_C>=L_N. It can also require connected spanning circulation or a specified cycle-rank threshold. Neutral rewiring preserves these variables, so [G,R]=0. A connected spanning circulation support requires L_N>=N. Fully circulating Q-edge graphs have e_C=Q.

The main input is the same coherent uniform neutral ray. Its conserved full-law moments are

\[
E=\langle H\rangle=\kappa Q,\qquad
\sigma^2=\|(H-E)\psi(T)\|^2=2Qt^2.
\]

The results below concern this energy/variance resource and closed evolution. They do not automatically apply to a different initial energy, postselection, driving, or an incoherent neutral mixture with a different rewiring energy.

## 2. Hopping norm resolved by vertex number

In rotor-angle space put c_i=sum_e B_ie(theta)b_e, with one unit-modulus coefficient per present half-leg, including the frozen sign. For each graph tr(BB^dagger)=2Q. Fermion exclusion gives sum_i c_i^dagger c_i<=2Q as an operator. This is pointwise in rotor angles and remains valid after imposing Gauss constraints.

On the n_v=m sector, the row creation map (v_i) -> sum_i a_i^dagger v_i has squared norm m+1, since its product with its adjoint on the output sector is sum_i n_i=(m+1)I. Therefore

\[
\boxed{\|A\Pi_m\|\le\sqrt{2Q(m+1)},\qquad
\|A^\dagger\Pi_m\|\le\sqrt{2Qm}.}
\]

The second follows by adjunction from the first. Orthogonality of different output vertex-number sectors and a triangle inequality give

\[
\boxed{\|K\Pi_{n_v\le k}\|\le
|t|\sqrt{2Q}(\sqrt{k+1}+\sqrt k)=b_k.}
\]

There is an exact improvement in vacuum. All edge modes are occupied, so a hop is necessarily edge-to-vertex. Its target identifies its parent uniquely: the sole hole edge and occupied vertex identify the half-leg whose electric shift is reversed. Distinct vacuum parents therefore have disjoint hopping targets, each with exactly 2Q contributions of magnitude |t|. Thus

\[
\boxed{\Pi_0 K^2\Pi_0=2Qt^2\Pi_0.}
\]

This identity is on the full unbounded physical space, not just a numerical projection. The direct checker enumerates all outgoing targets for 105 raw vacuum parents and verifies 840 distinct targets at N4,Q4.

## 3. A probability bound uniform in time

Let P=G Pi_(n_v<=k), and p(T)=||P psi(T)||^2. Since D and R commute with P and R>=0,

\[
P(D+R)P\ge\epsilon_E L_N P.
\]

Project the exact equation (H-E)psi=(D+R+K)psi and take norms. With delta=epsilon_E L_N,

\[
\delta\sqrt p\le\|P(H-E)\psi\|+\|PK\psi\|
\le\sigma+b_k.
\]

Consequently

\[
\boxed{\sup_T\Pr_T(G,\ n_v\le k)
\le\min\left(1,
\frac{2Qt^2}{\epsilon_E^2 L_N^2}
(1+\sqrt{k+1}+\sqrt k)^2\right).}
\]

For Q=nu N, L_N=sN and k=k_N=o(N), the bound is O((k_N+1)/N) and tends to zero. It holds even for times chosen as functions of N. **A nonvanishing extensive-circulation phase from this input cannot carry only a subextensive number of vertex fermions.** This is not a no-go theorem for an unrestricted charged phase.

At Q=2N with every Q edge circulating, frozen epsilon_E=1/4 and t=1 give

\[
\Pr_T(G,n_v\le k)\le
\min\left(1,\frac{16}{N}(1+\sqrt{k+1}+\sqrt k)^2\right).
\]

The constants are loose, especially for small N or k growing almost linearly. These are upper bounds, not predicted probabilities. The argument does not fit the decreasing N4-N6 populations or extrapolate their Taylor coefficients.

### Vacuum refinement

For P=G Pi_0, PKP=0 and ||PK(1-P)||<=b=|t|sqrt(2Q). Since PHP-E P>=delta P,

\[
\delta\sqrt p\le\sigma+b\sqrt{1-p}.
\]

For delta<=sigma use p<=1. Otherwise solving this inequality gives

\[
p\le\left[\frac{\delta\sigma+b\sqrt{\delta^2+b^2-\sigma^2}}
{\delta^2+b^2}\right]^2.
\]

Here sigma=b. For fully circulating Q=2N graphs this becomes

\[
\boxed{\sup_T\Pr_T(G,n_v=0)\le
\frac{64N}{(N+16)^2}\quad(N>16).}
\]

For the broader connected-spanning event, delta>=N/4 gives 256N/(N+64)^2 for N>64. The formulas are analytical rational functions of N. For example the fully circulating vacuum bound is about 0.060592 at N1024 and 0.000061034 at N1048576. It does not bound the sum of all charged and vacuum circulation sectors by the same expression.

## 4. Mean energy requires matter but does not prohibit the charged phase

Cauchy-Schwarz gives |<A>|<=sqrt[2Q<n_v>]. The electric-capacity inequality and R>=0 then imply, with x=<n_v>,

\[
\epsilon_E L_N\Pr(G)+(2m+\epsilon_E/2)x
\le2|t|\sqrt{2Qx}.
\]

For Q=nu N, L_N=sN and rho=x/N this is

\[
\boxed{\epsilon_E s\Pr(G)+a\rho
\le2|t|\sqrt{2\nu\rho},\qquad a=2m+\epsilon_E/2.}
\]

When the discriminant is nonnegative, probability at least p requires

\[
\rho\ge\left[\frac{|t|\sqrt{2\nu}
-\sqrt{2\nu t^2-a\epsilon_E s p}}a\right]^2.
\]

For a fully circulating Q=2N event of probability one, this coarse bound requires at least about 0.018125 mean vertex fermions per vertex. It is a necessary bound, not an optimum. It includes negative hopping and does not mistake positive electric energy for total energy.

## 5. A capacity-valid charged state on any supplied four-regular host

Let G_host be any connected simple four-regular graph, Q=2N. Orient an Euler tour; each vertex has two positive and two negative half-legs. Start with all edge fermions occupied, vertex vacuum, and fields (+1,-1) on each directed edge. The directed support is Eulerian and connected, so every edge lies on a directed cycle.

Assign each edge mode to its positive endpoint. At vertex i the two assigned edge modes and a_i have a three-state physical subspace: no hole/empty vertex, or one of two holes/occupied vertex. Relative to the unit-flow baseline, each hole changes its positive field from 1 to 2. Gauss's law holds because the vertex occupation equals the number of its assigned holes. The local vertex mode permits at most one hole in that pair.

Use the local state

\[
|\omega_i\rangle=\sqrt{1-\eta}|0\rangle
-e^{i\phi}\sqrt\eta|B_i\rangle,
\]

where |B_i> is the normalized bright superposition of the two positive-half-leg hopping targets, with the frozen fermion/orientation signs. The N even-parity local blocks commute. Their product is a normalized physical pure state with exactly M total fermions. Every configuration has r=0 and every graph edge still circulates. Its support has 3^N configurations, E^2 between 4N and 7N, and one full H action includes all outgoing fields, up to E^2=7N+3. No electric cutoff is introduced in the small checks.

The unit-flow subspace compression is a sum of independent local stars,

\[
h_i=\begin{pmatrix}0&\sqrt2t&0\\\sqrt2t&a_c&0\\0&0&a_c\end{pmatrix},
\qquad a_c=2m+3\epsilon_E/2=19/8.
\]

Define b_eta=sqrt[2eta(1-eta)] and

\[
\ell=a_c\eta-2t b_\eta\cos\phi,\qquad
q=2t^2+a_c^2\eta-2a_c t b_\eta\cos\phi.
\]

The full-H mean and variance are

\[
\boxed{\langle H\rangle/N=1+1/2+\ell,}
\]
\[
\boxed{\operatorname{Var}(H)/N=q-\ell^2
+2t^2(1-3\eta/2+\eta^2).}
\]

The last term is a full-law outgoing residual, not an omitted correction. Negative-half-leg hops leave the unit-flow band. Each target's anomalous negative field identifies the responsible edge, so outputs from different edges are orthogonal. The squared leakage per edge is t^2[(1-eta/2)(1-eta)+(eta/2)eta]. Local product structure gives the formula above. Rewiring annihilates every input configuration.

The formulas depend on degree and local state parameters, not host dimension. Optimizing the local star gives eta about 0.178474 and mean energy 0.840839N, below the neutral input's N. This is a variational comparator, not the exact full-H ground state. Choosing phi=0 and

\[
\eta_0=\frac{60-32\sqrt2}{291}\approx0.0506707
\]

gives exactly <H>=N. Its variance is N(41/16-3eta_0+2eta_0^2), about 2.415623N, so it does **not** match the neutral input's variance 4N. Rational arithmetic in Q(sqrt2) verifies the mean and local moments.

## 6. Even mean and variance matching do not select dimension

A further supplied pure-state construction matches both moments. Set eta=1/4. Choose opposite Euler orientations and phases with

\[
\cos\phi_\pm=\frac{35/32\mp h}{\sqrt{3/2}}.
\]

The two product branches have mean energies N(1+/-h). Their supports, and their one-H-action ranges, are orthogonal: reversing every edge changes signs on Q edges, whereas one hop changes just one half-leg. Thus an equal coherent superposition has

\[
\langle H\rangle=N,\qquad
\operatorname{Var}(H)=N[31/16+(N-1)h^2].
\]

Set h^2=33/[16(N-1)]. For every N>=128 both phase cosines lie in [-1,1], as checked by exact rational inequalities. This gives

\[
\boxed{\langle H\rangle=N,\quad\operatorname{Var}(H)=4N,
\quad\Pr(\text{full connected circulation})=1.}
\]

The state is pure, with <n_v>=N/4 and <E^2>=19N/4. Every connected simple four-regular host admits the construction. Supplied one-dimensional ring-power, two-dimensional square-torus and three-dimensional diamond families therefore admit identical values of these resources at common sizes N=16k^6 (N>=128 for the two-branch match). Their intrinsic balls differ.

These supplied states are not shown to have the same entire spectral measure as the neutral ray, belong to its same pure permutation sector, or be reachable from it. Fixed hosts break label symmetry; a permutation twirl preserves the moments and circulation but changes purity. No pure symmetry-sector matching is silently inferred. The exact conclusions are energetic admissibility and insufficiency of mean/variance/flux resources alone for dimension selection. The model can still distinguish topology through higher spectral information and dynamics.

## 7. Persistence and the next spectral gate

For the single product family, full-circulation loss has initial probability coefficient Q t^2(1-eta/2)(1-eta). One hop cannot erase both fields on any initially through-flowing edge, so loss of full protection has zero quadratic coefficient and begins no earlier than quartic order. These are short-lag controls, not lifetimes or a proof of persistent geometry. Requiring every single edge to keep exactly the same circulation property is stricter than preserving a usable coarse metric.

The next probability calculation should address the **charged** target sector, with vertex number of order N, extensive cycle rank and a stated metric event. Its overlap with the initial ray's occupied spectral directions matters. Let Pi_init project onto the invariant cyclic subspace generated by bounded spectral functions f(H)psi_0; using the exact initial permutation sector instead gives a weaker valid upper bound. For an energy interval I of half-width W, put s_I=||P_G 1_I(H)Pi_init||^2 and let tau be the conserved initial weight outside I. Then

\[
\Pr_T(G)\le\min(1,(\sqrt{s_I}+\sqrt\tau)^2),\qquad
\tau\le2Qt^2/W^2.
\]

A useful upper bound needs nontrivial control of s_I and the actual initial spectral tail; neither is supplied by this packet. There is a specific reason to avoid an unrestricted moment-only version. For a G containing the matched comparator, including the full-circulation event, that comparator's mean N and variance 4N imply ||P_G 1_I(H)||^2>=1-4N/W^2. Combining that unrestricted norm with only the neutral ray's Chebyshev tail upper bound 4N/W^2 makes the displayed probability estimate vacuous. Initial spectral directions, sector restrictions or stronger tail information are necessary to improve that route. This argument does not assume that the comparator meets an untested persistence or propagation criterion.

If appreciable charged-network weight is established, test persistence, intrinsic ball growth and pre-/post-source relational propagation in the same unfiltered state, with recording costs and neutral shortcuts accounted for.

## 8. Verification scope

`checks.py` verifies exact rational identities, the unregulated vacuum Gram identity, physical raising-block bounds in an independently built 6,075-state matrix, and 120 time-evolution cases of the general residual bound in 24 separate dense matrices. `prepared.py` compares two implementations of complete unbounded H action on six supplied states through N8, plus the opposite-orientation pair identity. No nonzero amplitude pruning or electric-field cutoff is used. Floating amplitudes and norms are not formal interval-roundoff certificates.

`hosts.py` performs exact integer graph/BFS checks on nine supplied hosts through N11664. Their quantum moments at large N follow from factorization, not a large Hamiltonian solve or an emergent-dimension fit. `verify.py` freshly reruns every checker and validates pinned sources and outputs. Independent external review of the new proofs remains pending.

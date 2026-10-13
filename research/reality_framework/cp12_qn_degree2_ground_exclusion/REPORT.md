# CP-12 Q=N: a ground-state exclusion for asymptotically degree-two geometry

**13 October 2026. Independent proof candidate, awaiting Claude Code/Work audit.** Additive research only: unchanged CP-12 Hamiltonian, unbounded integer U(1) electric rotors, exact Gauss constraints, all graph mixing, frozen (kappa,m,t,epsilon_E,Gamma)=(1/2,1,1,1/4,1/10). Neither GRUT-RAI nor Stage-3 Card 3 altered. Parent GRUT III adjudication: 694463f.

## Main result

In Q=N let D_2(G)=sum_i |degree_i(G)-2| and P_K include *all* physical matter/gauge states on graphs with D_2<=K. For **any K_N=o(N)** and any full-Hamiltonian normalized ground state psi_N, the argument below yields

\[
\boxed{\langle\psi_N|P_{K_N}|\psi_N\rangle\longrightarrow0.}
\]

For K=0 the quantitative estimate improves to O(N^-3); spanning rings and disjoint unions of cycles alike are excluded as dominant full ground-state geometry. Not an exclusion of extensive-defect geometries, non-equilibrium prepared rings, other Q, gravitational dynamics or all CP-12 spacetime phases.

## 1. A two-star, unregulated local quantum lower bound

The fixed-graph compression includes rewiring **diagonal** a_N=Gamma(M-Q)/M=(1/10)(1-2/(N-1)). Split the full fixed-graph block as H_G=sum_i h_i: each h_i gets vertex mass 2m n_i, all incident hopping, and half each incident edge's cost, electric energy and neutral-removable rewiring diagonal. This is the exact parent star decomposition in cp12_remainder_and_mixing/DERIVATION.md.

For any adjacent vertices i,j of degree exactly two, h_i+h_j involves only two central vertex fermions and three incident edge fermions. Relax Gauss constraints on outside vertices but keep both central vertex Gauss constraints and the three edge-center Gauss constraints. The resulting larger local space has five occupation bits and one free unbounded integer electric value z. Each of the four local incidence hoppings preserves z. The local matrix in every z sector is 32x32; no fixed loop holonomy or outer-vertex charge state is selected.

At a=19/200, an **exact rational LDL/Schur elimination** proves all 32 pivots of H_pair(z)+27/125 I positive for each integer -5<=z<=5. The smallest floating diagnostic pivot is about 0.00076436767; positivity is established by rationals, not eigenvalue floating-point comparisons. For |z|>=6, all six electric legs equal +/-z plus offsets at most 3, with sum of electric-energy weights =1/2. Thus electric energy >=(1/2)(|z|-3)^2>=4.5; edge cost=1, masses and removable-edge diagonal nonnegative, and the operator norm of four fermion–rotor hopping bilinears <=4. So these sectors have H_pair>=1.5. Consequently, with **unbounded** integer fields,

\[
\boxed{h_i+h_j\ge -27/125=-0.216\quad(a\ge19/200).}
\]

Monotonicity in a follows from the positive neutral-removable-edge projector. Since a_N>=19/200 for N>=41, this applies to all sufficiently large frozen-Q=N graph compressions. Even a triangle, whose two outer endpoints coincide, is covered: disregarding that exterior vertex's Gauss constraint only enlarges the local space.

For any exactly degree-two graph G, including arbitrary unions of differently sized cycles, every vertex is contained in two such adjacent pairs and therefore

\[
\boxed{\inf\sigma(H_G)\ge(-27/250)N=-0.108N.}
\]

## 2. Independent full-H quantum variational state

Use the already independently checked Gauss-invariant **coherent incidence-matching** trial of cp12_remainder_and_mixing/DERIVATION.md, section 4. It lies in the positive rewiring Laplacian's kernel, superposes all compatible graph configurations and has exact tridiagonal compressed Hamiltonian

\[
T_{kk}=\kappa Q+\tfrac{17}{8}k,\qquad
T_{k,k+1}=t(k+1)\sqrt{\frac{Q-k}{M-k}\frac{m_{k+1}}{m_k}}.
\]

Here m_k=binom(N,k)H_k, H_0=1, H_1=N-1 and H_k=(N-1)H_(k-1)-(k-1)H_(k-2). The derivation holds equally for Q=N; no tuning or external flux inserted.

An alternating-sign tent wavepacket centered at k=N/7 with width floor(sqrt(N)) has limiting full-H Rayleigh energy density

\[
u(1/7)=45/56-(12/7)\sqrt{2/7}
=-0.112752829414\ldots<-11/100=-0.11.
\]

The last strict comparison has an exact rational proof: 288/343-(45/56+11/100)^2=69113/13720000>0. Thus for sufficiently large N the **full** ground energy E0<-0.11N, separated extensively from every exactly degree-two compressed graph block. This is *not* a graph-by-graph perturbation theory.

The independent finite certificate using integer matching recurrences, rational squared couplings and outward Decimal evaluations gives:

| N | Coherent trial E/N upper | Degree-2 compression E/N lower | Ground exact-degree-2 probability upper |
|---:|---:|---:|---:|
| 42 | -0.1107852432 | -0.108 | 0.101724 |
| 64 | -0.1112796915 | -0.108 | 0.020931 |
| 128 | -0.1120794346 | -0.108 | 0.001706 |
| 256 | -0.1124767401 | -0.108 | 0.000177742 |
| 512 | -0.1125655797 | -0.108 | 0.000021405 |
| 1024 | -0.1126663971 | -0.108 | 0.000002564 |

These are bounds from small tridiagonal trial matrices, **not** exact full-H diagonalizations at these N.

Let P_2 include every exactly degree-two graph in Q=N. No single edge exchange connects two different 2-regular graphs. The independently established scalar rewire-cut bound, which dominates every Gauss-constrained physical fiber, is

\[
b_N^2=\|P_2H(1-P_2)\|^2\le\frac{\Gamma^2}{M^2}3Q(M-Q)=O(1/N).
\]

The full ground eigenequation and the extensive energy gap imply

\[
\boxed{\langle P_2\rangle_{\rm ground}\le b_N^2/(b_N^2+\delta_N^2)=O(N^{-3}).}
\]

## 3. Near degree-two extension, including internal graph mixing

For D_2(G)<=K, call a vertex with d!=2 defective; there are at most K and their summed degrees <=3K. Let E_good contain edges joining two degree-two vertices. The decomposition

\[
H_G=\tfrac12\sum_{ij\in E_{\rm good}}(h_i+h_j)+\sum_i c_i h_i
\]

has c_i=1-(# degree-two neighbors)/2 for degree-two vertices, and c_i=1 at defective vertices. Every star has h_i>=-|t|sqrt(d_i), by positivity of its diagonal and the normalized central fermion hopping anticommutator. The residual sum therefore has lower bound -[(3sqrt2)/2+4]K >=-7K. Since |E_good|<=Q=N,

\[
\inf\sigma(H_G)\ge-0.108N-7K.
\]

The full graph-rewiring offdiagonal inside the near-degree-two class is bounded by weighted Schur estimates. Any swap touching no defective vertex increases D_2 by 2 or 4. At most C_N(K)=3K(M-Q)+QK(N-1) swaps touch a defective vertex. With weights eta^{D_2(G)},

\[
\|P_K V P_K\|\le
r_N(K)=\frac{\Gamma}{M}\min\left\{Q(M-Q),
\inf_{0<\eta\le1}[\eta^2Q(M-Q)+\eta^{-4}C_N(K)]\right\}
=O(N^{2/3}K^{1/3})=o(N)
\]

for K=o(N), with r_N(0)=0. The crossing norm obeys

\[
\|P_KH(1-P_K)\|^2\le\frac{\Gamma^2}{M^2}Q(M-Q)C_N(K+4)
=O(N(K+4)).
\]

Every full-ground eigenstate obeys the sector-resolvent inequality with energy gap delta_N>=0.002N-o(N), yielding

\[
\boxed{\Pr_{\rm ground}(D_2\le K_N)=O((K_N+4)/N)\to0\quad(K_N=o(N)).}
\]

## 4. Independent quantitative controls and scope

An independently coded cycle Gauss/holonomy Hamiltonian reproduces the owner's published N=6 fixed-graph block energy -0.635203537847... at a_6=.06. For N=9 and a_9=.075, it finds E(C9)-3E(C3)=+0.003376757544..., with agreement across electric cutoffs three and four. These are corroborating small-system calculations, not large-N proof. A separate exhaustive N=6 graph-switch control checks the near-degree-two weighted counting.

**Reproducibility packet:** \`GRUT_CP12_QN_Degree2_Ground_Exclusion_2026-10-13.zip\` contains the full report, \`exact_pair_certificate.py\`, \`cluster_bound.py\`, \`independent_cycle.py\`, \`outward_ground_certificates.py\`, \`matching_trial.py\`, \`near_degree2_count_check.py\`, recorded JSON and SHA256 manifest. All listed checks ran in the independent container. This GitHub report is a transparent mathematical handoff; partner should execute the packet verbatim.

**Status:** candidate theorem with exact rational finite-star certificate, proven operator-norm/combinatorial steps, reproducible variational computations. Peer/partner review remains essential, especially the embedding of the local two-star matrix into *all* physical Gauss sectors and the finite-sector Schur estimates. If confirmed, retire Q=N *ground-state near-degree-two geometry* as a CP-12 spacetime route, not the entire frozen Hamiltonian.

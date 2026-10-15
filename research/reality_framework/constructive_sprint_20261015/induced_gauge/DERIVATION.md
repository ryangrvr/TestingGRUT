# Virtual matter generates gauge-loop dynamics in CP-12

This calculation adds no term to CP-12. It derives a Wilson-loop coupling from
the existing signed matter hopping. The exact frozen-parameter calculation is
the complete three-vertex, three-edge graph; there are no absent edges, so the
frozen rewiring operator is identically zero in this benchmark. This is a
mechanism check, not a graph-selection experiment.

## 1. Projection identity, energy conventions, and physical scope

Write

\[
\widetilde H=H-\kappa Q=2m n_v+F+R+K,\qquad
F=\frac{\epsilon_E}{2}\sum_{ie}E_{ie}^{2}.
\]

Let \(P=1_{n_v=0}\), \(\mathsf Q=1-P\), and let every occurrence of the scalar
\(z\) below refer to the energy of \(\widetilde H\), **excluding** the constant
\(\kappa Q\). It does not denote the energy of the full Hamiltonian. Whenever
\(z-\mathsf Q\widetilde H\mathsf Q\) is invertible,

\[
H_{\rm eff}(z)=P\widetilde HP+
P\widetilde H\mathsf Q
(z-\mathsf Q\widetilde H\mathsf Q)^{-1}
\mathsf Q\widetilde HP,
\]

\[
P(z-\widetilde H)^{-1}P=(z-H_{\rm eff}(z))^{-1}.
\]

These are exact Feshbach projection identities, not a truncated
Schrieffer-Wolff series. The energy dependence encodes the eliminated sector's
memory. At the frozen mass and hopping, replacing this by an autonomous,
energy-independent vacuum Hamiltonian is an additional approximation that is
not made here.

The projection method is established physics: see Feshbach,
*Annals of Physics* **5**, 357–390 (1958),
[doi:10.1016/0003-4916(58)90007-1](https://doi.org/10.1016/0003-4916(58)90007-1).
For the controlled perturbative framework, see Bravyi, DiVincenzo and Loss,
[arXiv:1105.0675](https://arxiv.org/abs/1105.0675). The CP-12 coefficients below
are calculations of this packet, not results attributed to those papers.

## 2. Exact minimum-order triangle coefficient

Let \(|0\rangle\) be the triangle's vertex vacuum, with all three edge
fermions occupied and electric fields zero. Let \(W\) increase the electric
field by one around the oriented triangle: one half-leg increases by one and
the other decreases by one on each edge. It preserves all Gauss constraints.
The state \(W|0\rangle\) has six nonzero half-leg fields and returns to
\(n_v=0\).

Each hop changes one half-leg field by one. Six hops are necessary. Every
minimal path changes each of these six fields once. Before the sixth hop,
no proper subset of a simple cycle can return to the vertex vacuum while
obeying Gauss law. Thus every intermediate state lies in \(\mathsf Q\), exactly
as required for the resolvent expansion.

There are 48 allowed signed paths. Twelve have vertex populations
\((1,2,1,2,1,0)\); thirty-six have \((1,2,3,2,1,0)\).
The product of the physical hopping signs is negative on every path. The
unweighted full-Hamiltonian sixth-power matrix element is consequently

\[
\langle W0|K^6|0\rangle=-48t^6.
\]

For a geodesic path after hop \(h\), \(\sum E^2=h\). Set

\[
\begin{split}
A_1&=2m+\epsilon_E/2,&A_2&=4m+\epsilon_E,\\
A_{3a}&=2m+3\epsilon_E/2,&A_{3b}&=6m+3\epsilon_E/2,\\
A_4&=4m+2\epsilon_E,&A_5&=2m+5\epsilon_E/2.
\end{split}
\]

The coefficient in the exact resolvent's hopping expansion is

\[
\boxed{
\langle W0|H_{\rm eff}(z)|0\rangle
=t^6\,
\frac{12/(A_{3a}-z)+36/(A_{3b}-z)}
{(A_1-z)(A_2-z)(A_4-z)(A_5-z)}
+O(t^8).}
\]

The five negative denominators cancel the negative fermionic loop sign.
For \(m=1,\epsilon_E=1/4,z=0\),

\[
\boxed{\langle W0|H_{\rm eff}(0)|0\rangle
=\frac{65536}{653429}t^6+O(t^8).}
\]

The coefficient is exact rational arithmetic. It is the Taylor coefficient of
a parameter-family resolvent at small \(t\); substituting \(t=1\) into this
single term is **not** a controlled frozen-law approximation.

## 3. General loop-length law in the controlled static/heavy-matter regime

On a fixed graph, introduce gauge angles and the signed incidence matrix
\(B(U)\). Each edge column has two unit-modulus entries with opposite signs,
as prescribed by CP-12. With electric evolution and graph rewiring suppressed
for this static band calculation, the one-particle Hamiltonian is

\[
h(U)=\begin{pmatrix}2mI&tB(U)\\tB(U)^\dagger&0\end{pmatrix}.
\]

Filling the edge-derived low band with the prescribed number of fermions gives
the exact static energy function

\[
\mathcal E(U)=Nm-\operatorname{Tr}
\sqrt{m^2I+t^2L(U)},\qquad L(U)=B(U)B(U)^\dagger.
\]

Zero singular values contribute zero to this expression. A uniform convergent
series exists when \(t^2\|L\|<m^2\). Its beginning is

\[
\mathcal E=-\frac{t^2}{2m}\operatorname{Tr}L
+\frac{t^4}{8m^3}\operatorname{Tr}L^2
-\frac{t^6}{16m^5}\operatorname{Tr}L^3+\cdots.
\]

A simple oriented \(\ell\)-cycle contributes
\(\ell(-1)^\ell(W_\gamma+W_\gamma^\dagger)\) to
\(\operatorname{Tr}L^\ell\). Combining this with the binomial coefficient gives

\[
\boxed{
\mathcal E_{\rm magnetic}
=\sum_\gamma J_{|\gamma|}
(W_\gamma+W_\gamma^\dagger)+\cdots,\qquad
J_\ell=
\frac{\ell\,C_{\ell-1}\,t^{2\ell}}{(2m)^{2\ell-1}}>0,}
\]

where \(C_j=(2j)!/[j!(j+1)!]\) is a Catalan number. For example,
\(J_3=3t^6/(16m^5)\), \(J_4=5t^8/(32m^7)\),
\(J_5=35t^{10}/(256m^9)\), and
\(J_6=63t^{12}/(512m^{11})\).

Exact geodesic sums independently reproduce this law for \(\ell=3,4,5,6\).
Direct diagonalization of the static single-particle problem reproduces its
small-\(t\) flux dependence. The script also reports finite-electric diagonal
denominator controls for these loop lengths; for \(\ell>3\) those controls hold
the cycle graph fixed and omit rewiring. They must not be interpreted as full
frozen-CP-12 results at those larger sizes.

On a fixed finite graph, the displayed coefficient is also the leading
large-\(m\) virtual-loop coefficient when the electric and rewiring scales are
kept bounded compared with \(m\). Resolvent insertions of those terms supply
additional inverse-mass corrections. Uniformity at increasing graph size
requires a separate estimate. This packet does not establish that uniformity.

The sign is physically informative. For one independent loop, its term is
\(2J_\ell\cos\Phi_\gamma\), favoring holonomy angle \(\Phi_\gamma=\pi\).
Holonomy angle is conjugate to electric flux. A preference for \(\pi\)
holonomy is **not** a preference for a particular nonzero electric-flux
configuration, nor proof of persistent electric circulation. On overlapping
loops the holonomies are constrained and their preferences can compete.

At fourth order,
\(\operatorname{Tr}L^2=\sum_i d_i^2+2Q\): within this static, controlled limit,
virtual matter penalizes degree inhomogeneity at fixed edge count. That is a
concrete organizing contribution, but it selects neither a lattice nor three
dimensions. Many regular graphs have the same fourth-order energy.

### Executed static coordination test: every six-vertex, twelve-edge graph

The script evaluates all \(\binom{15}{12}=455\) labelled graphs at zero
incidence angles, \(t=1\), with two explicitly labelled heavy-matter controls,
\(m=4\) and \(8\). These values test the expansion; they do not replace the
frozen \(m=1\) law. For every graph \(\|L\|\le6\), so the expansion ratios
are at most \(.375\) and \(.09375\), respectively. Define

\[
\mathcal E_4=-\frac{Q}{m}
+\frac{\sum_i d_i^2+2Q}{8m^3},\qquad
\mathcal E_6=\mathcal E_4-\frac{\operatorname{Tr}L^3}{16m^5}.
\]

The alternating square-root series implies
\(\mathcal E_6\le\mathcal E\le\mathcal E_4\) in this regime. Every numerical
case passes this analytic bracket. The exact static band minima occur on all
15 labelled four-regular graphs, each the complement of a perfect matching
(the octahedral graph). The next degree class has
\(\sum_i d_i^2=98\), compared with 96 for the regular graphs.

| Static control | Regular minimum \(\mathcal E\) | Lowest nonregular \(\mathcal E\) | Maximum \(|\mathcal E-\mathcal E_6|\) |
|---|---:|---:|---:|
| \(m=4,t=1\) | \(-2.797239384646\) | \(-2.794442113143\) | \(0.007707763136\) |
| \(m=8,t=1\) | \(-1.471834284387\) | \(-1.471388416707\) | \(0.000070900223\) |

The regular graph has spectrum \(L:\{0,4,4,4,6,6\}\), giving the independent
closed-form band energy
\(5m-3\sqrt{m^2+4}-2\sqrt{m^2+6}\).
This supports an earned preference for uniform coordination in the stated
controlled static test; it is not an inserted degree penalty.

An angle is an argument of this static energy function. A simultaneous sharp
rotor-angle configuration is not a normalizable finite-electric-energy
preparation. No such preparation is used or claimed. Nor is a graph dynamics
\(R+\mathcal E(g,U)\) executed here: a genuine Born-Oppenheimer reduction of
the full graph problem would need graph-offdiagonal overlaps of the dressed
fermion states, the electric kinetic term, and induced Berry/derivative
couplings. None can be silently omitted at the frozen parameters.

## 4. Frozen-parameter full Hamiltonian and exact Feshbach check

The actual frozen parameters are \(m=t=1,\epsilon_E=.25,\kappa=.5,\Gamma=.1\).
The benchmark is \(N=3,Q=3\), where the graph is unique and \(R=0\) identically.
The physical Gauss-law basis is restricted only by the stated regulator
\(|E_{ie}|\le L\). Results were compared at
\(L=1,2,3,4,6,8,12,16\); at \(L=16\), the full basis has 638 states and
the vacuum has 33 flux states.

The limiting displayed values are

| Quantity | Frozen triangle result |
|---|---:|
| Ground energy excluding \(\kappa Q=1.5\) | \(-1.96729955318915\) |
| Full ground energy | \(-0.46729955318915\) |
| Lowest finite-regulator eigenvalue of \(\mathsf Q\widetilde H\mathsf Q\) | \(0.1381571415968\) |
| Ground-to-eliminated-block spectral distance | \(2.105456694786\) |
| Ground vacuum weight | \(0.563578417645\) |
| Ground \(\langle n_v\rangle\) | \(0.508453854570\) |
| \(\langle W0|H_{\rm eff}(z_{\rm gs})|0\rangle\) | \(+0.0194943882733162\) |
| \(\langle W0|H_{\rm eff}(0)|0\rangle\) | \(+2.84590113171424\) |
| Ground \(\operatorname{Re}\langle W\rangle\) | \(-0.0456234862668\) |

The exact Feshbach eigen-equation residual is approximately \(2.1\times10^{-14}\).
The ground-energy change from \(L=12\) to \(16\) is about
\(3.1\times10^{-15}\), in ordinary float64 arithmetic. The positive effective
matrix element and negative Wilson expectation agree with the induced sign.
The ground state is not the evolved neutral input, and its substantial matter
admixture prevents treating the vacuum as a nearly closed physical system at
the frozen parameters.

The large difference between \(H_{\rm eff}(0)\) and
\(H_{\rm eff}(z_{\rm gs})\) is not a discrepancy: the resolvent is strongly
energy dependent. The eliminated block lies close to zero. Indeed, the single
sixth-order value \(0.1002955\) does not accurately approximate either
frozen-parameter matrix element. In the explicitly labelled hopping-family
check \(t=.5,.25,.125,.0625\), its relative error at \(z=0\) decreases through
\(.536,.104,.0247,.00609\), consistent with the next correction being of
relative order \(t^2\).

## 5. Control of the unbounded electric space

For the triangle, Fourier transformation of the electric shifts writes the
hopping operator as a direct integral of fermionic quadratic Hamiltonians.
Its norm obeys

\[
\|K\|\le\sum_{a=1}^{3}s_a(B)
\le\sqrt{3\operatorname{Tr}BB^\dagger}=\sqrt{18}.
\]

Outside the \(L=16\) field cube, at least one half-leg has magnitude at least
17. Therefore the omitted full space has
\(\widetilde H\ge 17^2/8-\sqrt{18}=31.8823593\), and the omitted matter
space has \(\mathsf Q\widetilde H\mathsf Q\ge33.8823593\).
Explicit integer enumeration of outgoing hopping destinations gives
\(\|B_{\rm boundary}\|^2\le\|B\|_1\|B\|_\infty=3\).
There are no omitted rewiring or diagonal matrix elements.

If a retained block has lower bound \(a\), its omitted partner lower bound
\(d>a\), and the cross-block norm squared is at most \(b^2\), then

\[
H\ge\frac{a+d-\sqrt{(d-a)^2+4b^2}}2.
\]

Using the eliminated-matter block gives

\[
\boxed{\mathsf Q\widetilde H\mathsf Q\ \gtrsim\ 0.0494859853>0.}
\]

The inequality's analytical form is exact; the displayed value uses numerical
retained eigenvalues. This controls both resolvents at \(z=0\) and
\(z=z_{\rm gs}\) without assuming that the cutoff is the law.

For a computed retained resolvent vector
\(x=(z-H_{\mathsf QQ,L})^{-1}q_0\), its full-space residual is the retained
linear-solve residual together with \(B_{\rm boundary}x\). If the full eliminated
block has lower bound \(a_*\), then

\[
|\langle q_1,(z-H_{\mathsf QQ})^{-1}q_0\rangle
-\langle q_1,x\rangle|
\le\frac{\|q_1\|}{a_*-z}
\big(\|r_{\rm finite}\|+\|B_{\rm boundary}x\|\big).
\]

Here \(\|q_1\|=\sqrt6\). Evaluating this bound gives about
\(3.3\times10^{-13}\) at \(z=0\), and \(4.7\times10^{-16}\) at the computed
ground energy. The ground itself is isolated below \(-1.5\) by a min-max block
comparison: use
\(H\ge(A-\beta I)\oplus(d-b^2/\beta)\), choose
\(\beta=b^2/(d+1.5)\), and check that the retained second eigenvalue minus
\(\beta\) exceeds \(-1.5\). A Temple residual estimate then bounds the ground
energy using this gap.

All floating calculations and residual evaluations use ordinary float64.
These are analytically derived error controls with numerical evaluation,
not formal interval-arithmetic certificates. In particular, apparent tiny
Temple widths must not be advertised as certified digits beyond floating
roundoff.

## 6. Constructive implication and remaining task

CP-12 does contain an internal process that moves between Gauss-invariant
vacuum electric-loop configurations. It is generated by the physical
matter-hopping term; no magnetic plaquette term was inserted. The low-energy
description naturally contains Wilson loops and an explicit loop-length/sign
relation. This is a more structured mechanism than arbitrary independent edge
protection.

It does not select the supplied graph, generate three-dimensional volume
growth, or establish a causal cone. Rewiring remains all-to-all in the full
theory. The earlier all-time low-matter obstruction is unchanged: deriving an
effective vacuum gauge theory cannot establish dominant extensive circulation
reached from the specified neutral initial ray. Extensive charged networks
and their actual dynamical accessibility remain necessary research targets.

The next useful calculation is a **charged** low-energy/cyclic-sector
effective theory retaining the vertex matter, with the induced loop operator
and rewiring both kept. Its first test should compare spectral weights and
spatially resolved propagation on the same internally formed graph sector,
rather than equate this triangle benchmark with emergent space.

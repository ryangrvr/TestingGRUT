# Exact geometric reconstruction boundary

## Response of the spring Hamiltonian

Let $x,p$ denote the canonical displacement and momentum of $n$
unit-mass sites, and let $L$ be the positive weighted spring Laplacian.
A common pin gives $V=I+L>0$, with closed quantum Hamiltonian
$$
H=\tfrac12(p^{\mathsf T}p+x^{\mathsf T}Vx).
$$
The source is $-f(t)^{\mathsf T}Px$ for a coordinate selector $P$.
The Heisenberg equations for the driven system give
$\ddot x+Vx=P^{\mathsf T}f$. The retarded displacement susceptibility
is independent of the state because the canonical commutator is scalar:
$$
\begin{aligned}
\chi_P(t)
 &= \Theta(t)P V^{-1/2}\sin(tV^{1/2})P^{\mathsf T},\\
\chi_P^{(2j+1)}(0^+)
 &=(-1)^j P V^j P^{\mathsf T},\\
\widehat\chi_P(s)
 &=P(V+s^2 I)^{-1}P^{\mathsf T},\qquad \Re s>0 .
\end{aligned}
$$
The derivative is taken on the positive-time branch. The last line is
a Laplace transform. The calculator uses the positive real parameter
$z=s^2$ in its resolvent checks. The full projected spectral measure
has atoms at stiffness eigenvalues and weights
$P\Pi_\lambda P^{\mathsf T}$, with $\Pi_\lambda$ the orthogonal spectral
projector. This measure gives the susceptibility and its moments.

For fixed observed labels $a\ne b$, the first nonzero Taylor
derivative has order $2d(a,b)+1$, where $d$ is the unweighted hop
distance in the physical spring graph. For smaller powers there is no
walk between the sites. At power $d$ every contributing walk is a
shortest path and must use $d$ off-diagonal entries, all negative.
Consequently the derivative is the positive sum of products of spring
weights on the shortest paths. Thus the data determine the graph
distance **between measured sites** and all weights on direct
links between them. Complete calibrated site access gives $V$ from
the third derivative.

## Classification with few inaccessible sites

Let $W$ be the stiffness of another connected network with the same
pin, port labels and full Krylov span. After aligning measured labels,
equal responses mean $P V^j P^{\mathsf T}=P W^j P^{\mathsf T}$ for all
$j\ge0$. Define
$$
Q\, V^j P^{\mathsf T}u=W^j P^{\mathsf T}u .
$$
The inner product between the generators with indices $j,k$ and
vectors $u,v$ is $u^{\mathsf T}P V^{j+k}P^{\mathsf T}v$. Moment
equality gives equal Gram forms, so every linear relation between
the generators on one side also holds on the other. The map extends
to an orthogonal isomorphism of the two full site spaces. It fixes
the measured coordinate vectors at $j=0$ and satisfies
$QV=WQ$. Conversely every such similarity gives identical port moments
and susceptibility. In finite dimension, full access to the modes
is equivalent to the absence of a nonzero eigenvector vanishing at
every port. The orthogonal complement of the Krylov span is
$V$-invariant, so a nontrivial complement contains such an eigenvector.

Connectedness makes the eigenvalue 1 of $V$ simple, with eigenvector
$\mathbf 1$. Hence $Q\mathbf 1=\pm\mathbf 1$. A measured coordinate
fixed by $Q$ forces the positive sign. Order the coordinates with the
observed ones first. The remaining action is an orthogonal matrix on
the hidden block fixing its own constant vector. With zero or one
unmeasured sites this action is trivial. With two, its remaining
one-dimensional action is a sign, giving identity or the swap.
With at most two hidden sites the weighted graph is therefore
unique up to their relabeling. The conclusion uses unit masses,
common pinning, connectedness, exact data and full visibility.

## Exact nonisomorphic pair with three hidden sites

Measure site 0 in the following two graphs. The weights are in the
frozen units.

| Graph | Positive springs |
| --- | --- |
| Cycle | $01:27,\ 02:18,\ 13:2,\ 23:8$ |
| Pendant triangle | $01:12,\ 02:3,\ 03:30,\ 12:10$ |

The resulting Laplacians are
$$
L_{\rm cycle}=
\begin{pmatrix}
45&-27&-18&0\\
-27&29&0&-2\\
-18&0&26&-8\\
0&-2&-8&10
\end{pmatrix},
\qquad
L_{\rm triangle}=
\begin{pmatrix}
45&-12&-3&-30\\
-12&22&-10&0\\
-3&-10&13&0\\
-30&0&0&30
\end{pmatrix}.
$$
Their unweighted supports have zero and one triangle respectively.
They have the same edge count, total spring strength 55, masses,
pinning and normal-mode frequencies. For
$$
O=\tfrac13\begin{pmatrix}2&-1&2\\-1&2&2\\2&2&-1\end{pmatrix},
\qquad Q=1\oplus O,
$$
direct multiplication gives $O^{\mathsf T}O=I$, $O\mathbf 1=\mathbf
1$ and $Q L_{\rm cycle} Q^{\mathsf T}=L_{\rm triangle}$. The common
port is fixed. The matrix with columns $e_0,V e_0,V^2e_0,V^3e_0$ has
determinant of magnitude 1,807,920 on both graphs. Replacing $V$ by
$L=V-I$ changes these columns by a triangular transformation of
unit determinant. In particular the stiffness is cyclic from the
port, and the ambiguity cannot come from a mode added outside the
measured realization.

For any number of ports, append a common path of new, measured sites
from site 0. The transformation fixes all these coordinates and all
appended spring terms. A dark eigenvector of the extension would
vanish on the entire appended path and site 0, leaving a dark vector
of the original network. Full visibility and the distinct triangle
counts persist. If more than three hidden sites are wanted, attach
new hidden leaves to site 0 with distinct weights in $(100,101)$.
On the hidden block their diagonal eigenvalues are disjoint from
the core's hidden block, whose norm is less than 100. The scalar
coupling to site 0 is nonzero on every new leaf eigenvector and on
every core hidden-block eigenvector by the original rank result.
Thus no eigenvector can vanish on all ports. Three unmeasured
coordinates are the sharp threshold for this class, including as
the fraction of measured sites tends to one.

## Full port experiments and a separating probe

The orthogonal map of coordinates is a canonical transformation
acting only on the hidden degrees of freedom. It is implemented
on $L^2(\mathbb R^n)$ by $\psi(x)\mapsto\psi(Q^{\mathsf T}x)$.
This intertwines the two Hamiltonians and fixes the measured
position and momentum operators. Corresponding preparations,
including Gibbs states at every common inverse temperature and the
ground states, give the same projected correlations. For example,
$$
\langle xx^{\mathsf T}\rangle_\beta
=\tfrac12 V^{-1/2}\coth(\tfrac12\beta V^{1/2})
$$
in the frozen units. The same intertwining applies to all sequences
of controls and instruments on the fixed measured algebra. Equal
priors on the two networks then give maximum decision success $1/2$
from this access. Individually calibrated hidden sites fix a
different physical observable dictionary, so the graphs need
not agree on those new experiments.

The registered additional probe measures a force at site 3 through
the response of site 0. Its third derivative is $-V_{03}$, equal
to $0$ on the cycle and $30$ on the pendant triangle. The two
matrices have norm at most 91 by the absolute row sums. For a
time $t>0$ with $91t^2<42$, the Taylor remainder per network
after the cubic term is at most
$$
\frac{91^2 t^5}{5!}\frac{1}{1-91t^2/42}.
$$
The ratios of successive absolute terms starting at fifth order
are bounded by $91t^2/42$. At $t=1/20$ the signed response gap
is therefore greater than $29/50000$ in calibrated units.
A measured gap outside the certified interval would reject the
specified pair model or measurement calibration. The distinction
is predicted by the same ordinary spring law in both graphs.
[Exact certificate](RESULTS.json) · [Frozen input](WITNESSES.json)

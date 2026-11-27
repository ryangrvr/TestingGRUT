# Complementary interventions on the F1 springs

## Physical interface

Retain the F1 stiffness $V=I+L$, with its exact edges and unit
baseline masses. For a calibrated diagonal mass matrix $M>0$
and pin increment $D\ge0$, take the displaced Hamiltonian
$$
H=\tfrac12(p^{\mathsf T}M^{-1}p+x^{\mathsf T}(V+D)x)-f(t)x_0.
$$
At the same original port the retarded driving-point transfer at
positive Laplace variable $s$ is
$$
r_{M,D}(z)=e_0^{\mathsf T}(V+D+zM)^{-1}e_0,\qquad z=s^2.
$$
A static load setting is established before applying the force.
Independent physical tags name which oscillator couples to a
local control, without specifying its neighbors, distances,
or spring strengths. Additive pinning and inertial loading are
specified operators. The experiment presumes that either can be
changed without uncontrolled changes in the other terms. This
is a calibration resource, and the proof prices it in accessible
settings and supported coordinates. It presumes no hidden
displacement or force detector. [F1 witness](../grut_f1_operational_reconstruction_20261126/WITNESSES.json)

## What remains invariant and what separates the known pair

The F1 canonical change of variables fixes the port algebra.
With consistently transformed initial states, all port-only
control sequences remain equivalent, including nonlinear drive
and nonlinear readout there. An extra uniform mass multiplier,
a root-only pin shift, and any perturbation invariant under that
same transformation preserve the equivalence. In particular the
hidden block has nonzero off-diagonals throughout the F1
intertwiner, so a diagonal mass pattern is invariant under it
exactly when its three hidden entries are equal. An invariant
radial hidden interaction such as
$\gamma(\sum_{i=1}^3 x_i^2)^2/4$ also preserves equivalence.

For $G=(V+zI)^{-1}$, a calibrated pin $d>0$ at tag $i$ gives
$$
 r_{\varnothing}-r_i(d)=\frac{d G_{0i}^2}{1+dG_{ii}}.
$$
An additional mass $\mu>0$ at that site gives the same formula
with $d=z\mu$. Thus these distinct physical load channels
coincide at one shift if the effective impedances coincide.
The registered unit pin and unit mass at $z=1$ both give
$$
r_{3,\mathrm{cycle}}=\frac{2297}{17077},
\qquad
r_{3,\mathrm{triangle}}=\frac{715}{5591}.
$$
Their exact corresponding contrast is $632472/95477507$.
When the hypotheses are the frozen pair, one additional calibrated
intervention and one scalar value suffice. With zero such
interventions all port-only possibilities remain equivalent,
so this is minimal in the declared pair test.

The diagnosis is available even with the hidden tag identity
unknown within the pair. At large $z$, $G_{0i}$ starts at order
$z^{-2}$ with coefficient $w_{0i}$ for a neighbor of the port.
For cycle tag 3 the leading order is $z^{-3}$ with coefficient
198. The positive neighbor strengths are 27 and 18 in the cycle
and 12, 3, and 30 in the triangle. They have disjoint
positive squares. With a calibrated positive load, their
root-response losses have disjoint leading powers or
coefficients across the two hypotheses. At the two registered
finite shifts, all nine cross-hypothesis tag choices are
already distinct for each load type. This statement is
about the given pair. If the load strength is free at one
shift, both loss curves cover a common interval starting
at zero, since their ranges are $(0,G_{0i}^2/G_{ii})$ for
arbitrary positive $d$. A calibrated strength or richer
data are essential to the single-shift comparison.
[Exact values](RESULTS.json)

A finite-time version of the $z=1$ comparison truncates the
Laplace measurement at $T=20$. In the two unit-mass or
unit-increment settings, $I\le M\le 2I$ and $V+D\ge I$.
Thus the impulse response norm is at most $\sqrt{2}<2$.
At $s=1$ the two-tail error is at most $4e^{-20}$, bounded
by $4(2/5)^{20}<4.40\times10^{-8}$. The separation remains
greater than $0.0066242$ in the frozen units before other errors.

## Full reconstruction using paired loads

**Theorem.** For any connected positive spring network with
specified equal baseline mass and tag labels, a single scalar
$G_{00}$ plus two positive singleton loads per hidden tag and
one simultaneous positive load per unordered pair determine the
whole stiffness at a single fixed $z>0$. The baseline system
can have dark modes. The measurements number
$1+2(n-1)+\binom{n-1}{2}$. This is an explicit sufficient
setting count. The lower bound below concerns simultaneous support
size.

Write $u_i=G_{0i}>0$ and $v_i=G_{ii}>0$. Strict positivity of
$G$ follows from the nonnegative walk expansion of the inverse
of a strictly diagonally dominant spring matrix and connectivity.
For instance with $K=T-A$, where $T$ is its diagonal and
$A\ge0$ contains the springs, $\|T^{-1}A\|_\infty<1$ and
$K^{-1}=\sum_{k\ge0}(T^{-1}A)^k T^{-1}$. This supplies the
sign of $u_i$. The function
$$
\frac{d}{r_\varnothing-r_i(d)}=\frac{1+d v_i}{u_i^2}
$$
has known slope and intercept from the two measurements. Recover
$u_i$ from the positive square root and $v_i$ from the slope.

For a simultaneous pair $i,j$ let $\alpha=v_i+1/d_i$,
$\beta=v_j+1/d_j$ and $g=G_{ij}$. Woodbury reduction gives
the loss
$$
F_{ij}(g)=\frac{\beta u_i^2+\alpha u_j^2-2g u_i u_j}{\alpha\beta-g^2}.
$$
A principal submatrix of the spring matrix has a nonnegative
inverse by the same walk expansion. The Schur identity for
removing coordinate $i$ gives $u_j-u_i g/v_i\ge0$.
Removing $j$ gives the other inequality. The physical interval
is therefore
$$
0\le g\le \min(v_i u_j/u_i,\ v_j u_i/u_j).
$$
This upper bound is at most $\sqrt{v_i v_j}<\sqrt{\alpha\beta}$.
On the entire interval,
$$
F_{ij}'(g)=
-\frac{2(\beta u_i-g u_j)(\alpha u_j-g u_i)}{(\alpha\beta-g^2)^2}<0 .
$$
Thus a single pair value selects a unique physical cross entry.
The magnitude of the slope is at least
$2u_i u_j/(d_i d_j\,\alpha^2\beta^2)$ conditional on the
exact singleton data. Consequently injectivity is compatible
with poor conditioning for weakly coupled tags. In the
infinitesimal version,
$\partial_{d_i}r=-u_i^2$ and
$\partial_{d_i}\partial_{d_j}r=2u_i u_jG_{ij}$ at zero
load, from differentiating the inverse. The finite theorem
requires no derivative measurements.

The recovered full positive matrix gives $V=G^{-1}-zI$.
The exact fixture solves the pair quadratic and selects the
root in the physical interval without rounding. If the
baseline root is cyclic, all-frequency singleton data alone can
already suffice, since they determine the full first column
$G(z)e_0$ and its moments, which span the whole space in that
case. The worst-case assertion below drops that visibility
premise.

## All singleton settings can fail in the larger class

Attach one port by unit springs to all 16 vertices of either the
rook or Shrikhande graph, with unit pin and mass everywhere.
Label hidden sites by $\mathbb Z_4^2$. Rook adjacency means a
common row or column for distinct sites. Shrikhande adjacency
uses increments $(\pm1,0),(0,\pm1),(1,1),(-1,-1)$.
Both hidden adjacency matrices satisfy
$A\mathbf 1=6\mathbf 1$ and $A^2=4I+2J$ with
$J=\mathbf 1\mathbf 1^{\mathsf T}$. One has a hidden
four-clique and the other has none. The root is the unique
degree-16 vertex, and both networks have 64 unit springs.

The hidden block of $V+zI$ is $B=(z+8)I-A$ and the root
off-diagonal is $-\mathbf 1^{\mathsf T}$. With $w=z+8$,
$$
B^{-1}=\frac{wI+A}{w^2-4}
       +\frac{2J}{(w-6)(w^2-4)} .
$$
Both have the same diagonal and row sum $1/(z+2)$.
Hence block inversion gives
$$
G_{00}=\frac{1}{z+17-16/(z+2)},\qquad
G_{0i}=\frac{G_{00}}{z+2},\qquad
G_{ii}=(B^{-1})_{ii}+\frac{G_{00}}{(z+2)^2}.
$$
The entire singleton loss formula is identical for each
tag, all positive shifts and all nonnegative singleton
strengths, including singleton mass plus pin loads. A further
root-only load maps the common scalar $r$ to $r/(1+d_0 r)$,
so it preserves equality too. The off-diagonal Green entry differs for
$(0,0),(0,2)$ since the coefficient of $A$ in the displayed
inverse is positive. The double-load comparison gives
an exact gap $77/5742480$ at $z=1$. Because the singleton
inputs to $F_{ij}$ are identical and the cross entry differs,
strict monotonicity entails the gap for each positive $z$
and positive pair load. With this matched control,
simultaneous support size two is both sufficient and necessary
in the worst case of the stated static loading class.
[Registered replay](RESULTS.json)

## Nonlinear and independent interaction controls

Attach an independently observed oscillator $y$ with the
specified elastic term $(y-b^{\mathsf T}x)^2/2$ and its own
unit mass and pin. A collective vector $b=(0,1,1,1)$ is fixed
by the F1 change of variables. Its entire extended two-port
transfer is invariant. The independently calibrated local
$b=e_3$ is separating in the exact replay at both shifts.
An extra physical channel can share the old symmetry even when
it has genuinely new dynamics.

For the nonlinear test, start from corresponding harmonic
ground states, then change the Hamiltonian to
$H_\gamma=H_0+\gamma x_i^4/4$ with a root momentum kick $a$.
Use the coefficient linear in $\gamma$ and cubic in $a$ in the
expectation of $x_0(t)$. The kick shifts the harmonic solution
for site $i$ by $a K_{i0}(t)$, with
$K(t)=V^{-1/2}\sin(t V^{1/2})$. The first variation from the
quartic force is a retarded integral of the harmonic cubic.
Equivalently the quartic interaction commutator with a harmonic
displacement is a cubic times a scalar. These are identities
for the perturbative coefficient on the harmonic preparation,
where polynomial moments are finite. The term cubic in the
kick is state independent:
$$
C_i(t)=-\int_0^t K_{0i}(t-u)K_{i0}(u)^3\,du.
$$
For a leading propagation term $c t^p/p!$, the first coefficient
is $-c^4(3p)!/[(p!)^3(4p+1)!]$ at time power $4p+1$.
Tag 3 gives $(c,p)=(198,5)$ for the cycle and $(30,3)$ for the
triangle. Thus the leading cubic responses are
$-(131769/5788160000)t^{21}$ and $-(125/572)t^{13}$.
A calibrated tagged nonlinearity distinguishes them. Arbitrary
nonlinear port-only protocols still obey the earlier equivalence.
This calculation concerns the first variation, with no performance
claim for a finite large-amplitude nonlinear apparatus.

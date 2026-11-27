# GRUT-F2 independent internal review

A separately tasked internal AI built an exact replay after the freeze at
2db0de4c0c50aed3a7297964fc8ca9fedaa96991. It imports no production code.
The added access, dynamics and assumptions are checked against the frozen
protocol. External RAI and hardware replication remain open.
[Protocol](FROZEN_PROTOCOL.md) · [Parameters](PARAMETERS.json)

## Exact replay

The checker builds the shifted stiffness matrices directly and uses
rational $LDL^{\mathsf T}$ solves. It obtains both F1 networks' root
responses for every registered hidden unit pin and unit mass at $z=1,2$.
It also checks the root-only and uniform-mass controls, the two auxiliary
oscillator transfers, and the coefficient of the mixed quartic response.
The first invocation exposed integer division at the first LDL elimination
step. Forcing rational pivots removed that implementation error before the
registered replay was saved. The frozen inputs and gates did not change.

~~~bash
python3 critical_replay_f2.py --output /tmp/grut-f2-critic.json
~~~

The independent run took less than one second. Both networks were
reconstructed exactly from one baseline reading, two positive singleton
loads on every hidden tag and one positive loading of each hidden pair.
For the tag 3 contrast at $z=1$, the cycle gives $2297/17077$ and
the pendant triangle gives $715/5591$, in both the pin and mass tests.
The difference is $632472/95477507$. Across the three possible hidden
tags there are no cross-network collisions at either registered shift for
each perturbation type separately. This finite check gives no universal
latent-tag result for other graphs or unknown load strengths.
Read-only comparison with the production archive matched 274 common
scalars exactly, including the registered stiffness matrices and
the continuum signal floors.

For the time-truncated reading, the operator
$M^{-1/2}(V+D)M^{-1/2}$ has minimum eigenvalue at least $1/2$.
The root retains unit mass. Each impulse response therefore has amplitude
at most $\sqrt{2}<3/2$. From $e^{20}\ge20^{20}/20!$, the combined
missing Laplace tail at $s=1$ is bounded by the exact rational
$3\cdot20!/20^{20}=44547766263/640000000000000000$.
This weaker independent cap still certifies a gap greater than
0.0066242 in the registered units. The bound controls truncation only.
[Exact replay](critical_replay_f2.py)

## Rank-two inverse and the meaning of a tag

Positive springs and pinning give an irreducible positive definite
$M$-matrix at each shifted stiffness. Its inverse $G$ is strictly
positive from a nonnegative walk expansion. The positive root-to-tag
entry $u_i=G_{0i}$ is therefore recoverable from its square.
Two distinct calibrated positive singleton strengths recover $u_i^2$
and $d_i=G_{ii}$ from the rank-one loss formula.

For pair strengths $\lambda_i,\lambda_j$, set
$a_i=d_i+\lambda_i^{-1}$ and $a_j=d_j+\lambda_j^{-1}$. The loss as a
function of $c=G_{ij}$ has derivative

$$
-\frac{2(a_j u_i-cu_j)(a_i u_j-cu_i)}
       {(a_i a_j-c^2)^2}.
$$

Deleting tag $j$ leaves a spring principal subproblem with a
nonnegative inverse even when it disconnects. Its root-to-$i$
Schur formula gives $u_i-cu_j/d_j\ge0$. Deleting $i$ gives the
symmetric inequality. All physically possible $c$ lie on an interval
where the two derivative factors are strictly positive because the
loads are finite and $u_i,u_j>0$. The denominator is strictly
positive there as well. Hence the observed pair loss determines a
unique cross entry among admissible resolvents. The exact replay
selects that root of the pair quadratic. This recovers $G$ before
inverting it using calibrated load responses alone. The inverse can
be ill-conditioned under noise, particularly when a tag is weakly
coupled to the root.

The independently calibrated support and strength of each load supply
physical hidden coordinate axes. Enumerating the tags supplies their
count and a correspondence across settings. Uniform masses, known
pinning and time and readout calibration reduce the measured
$G^{-1}$ to the stated spring matrix. These inputs contain no spring
adjacency, but they are more access than F1's root-only interventions.
In particular a tunable actuator on a tagged oscillator is an
independent physical interaction.

## Singleton and nonlinear controls

The rook and Shrikhande graphs have parameters $(16,6,2,2)$
and hidden adjacency polynomial $(x-6)(x-2)^6(x+2)^9$.
The rook graph has eight four-cliques and the Shrikhande graph has
none. A uniform root is the unique vertex of degree 16 in either
rooted graph. Their common adjacency algebra gives
$A^2=4I+2\mathbf 1\mathbf 1^{\mathsf T}$. Block inversion then
fixes the same root, root-to-tag and tag diagonal Green entries
in both graphs. Each singleton loss depends only on those entries.
This gives exceptional all-frequency equality for every
nonnegative strength, independently of the finite replay. At
$z=1$ the pair setting gives $265/3408$ versus $131/1685$, a
gap of $77/5742480$. At $z=2$ it gives $9/136$ versus
$107/1617$. Identically calibrated root loads do not cure singleton
equality. A time sequence coupling one tag after another in a single
trial falls outside this static-reset obstruction. The sharp lower
bound concerns concurrent support size, with no claim of optimal
number of measurements.

The rooted pair defeats singleton-only reconstruction without a
visibility assumption. Full visibility would change the comparison.
If the root is cyclic, knowledge of the entire singleton functions
and calibrated tags determines the first inverse column and hence
the root-generated Krylov vectors, which identify the stiffness.
F2 retains the F1 result at the weaker original port and proves a
conditional reconstruction under independent perturbations.

The collective auxiliary oscillator is physically additional but
respects the F1 change of hidden coordinates. Its full two-port
response stays invariant. The auxiliary attached to site 3 changes
that response at both shifts. Likewise a root-supported nonlinear
drive or a radial hidden quartic respects the canonical equivalence
for corresponding preparations. With a tagged quartic and a root
kick, the cubic-in-kick part of the first variation is the convolution
of one retarded root-to-tag amplitude with its cube. The leading
powers are 21 and 13, with respective coefficients
$-131769/5788160000$ and $-125/572$ after factoring the quartic
strength and cube of the kick. In a Gaussian harmonic preparation
other kick powers can occur. The certified coefficients concern
the stated mixed derivative and do not by themselves certify a
finite nonlinear laboratory signal.

## Causal scalar-probe check

In the conformal variables the known mass and curvature couplings
expose the two potentials $X=\Omega^2-1$ and $Y=\Box\Omega/\Omega$
through $U_i=m_i^2 X+(1-6\xi_i)Y$. The conformal massless probe
has zero potential. Its entire restricted response agrees.
Equality through the early gate follows from a shortest
bump-to-detector arrival of $13/8$, after the early gate closes
at $3/2$. The future bump is outside the measurement's retarded
domain. These controls address causal support rather than
small-magnitude scoring.

The temporal and spatial supports of the light cones meet only
the source and late detector temporal plateaus on the full bump.
Their massless waves have scaled value $1/r$ by the shell theorem,
where $r$ is distance from the apparatus center. The massive
free tail satisfies $|J_1(x)|\le x/2$, giving the registered free
wave interval $[123/160,131/96]$. Integrating the null cone over
the bump's temporal width gives $1/32$, and the massive tail
contributes at most $m^2/6144$ in the selected conservative norm.
The perturbed source and advanced detector waves remain above
$3/4$ for the whole massive parameter interval. The resolvent
identity uses their product against the nonnegative potential
difference, so the late massive score strictly decreases in this
known bump family. The weak and strong magnitude floors are
$2.8610\times10^{-16}$ and $8.0466\times10^{-7}$.

For the minimal massless channel, one integration by parts moves
$\Box$ onto the product of the two free waves. On the bump
$\Box r^{-2}=2r^{-4}$, giving a negative linearized score
shift. The box derivative cap $32768$ follows separately from
the temporal and radial second derivative caps and the radial
term outside the interior plateau. The first Born remainder is
bounded by the slab norm and the deviation of $\epsilon\Box
b/(1+\epsilon b)$ from $\epsilon\Box b$. In the weak case its
ratio to the linear floor is below .07652, and the signed
magnitude lower bound is $4.6917\times10^{-15}$.
The strong minimal Neumann certificate is 256, leaving the
registered estimate inconclusive. This is distinct from an exact
null result.

Using the common nonnegative bump factor gives the independent
linearized ratio interval $[5.3642,46.187]$ for the magnitude of
minimal over massive response. This is an interval for the
specified channels and source profiles. If both effective
potentials were independently identified in space and time, the
pointwise rank-two relation and the consistency equation
$Y=\Box\sqrt{1+X}/\sqrt{1+X}$ would constrain the profile.
The present scalar scores do not identify those fields.
The specified $\xi$ is an independent calibration since
flat-region probes cannot recover it. Mass and gain calibration
also matters. For example a constant rescaling of $\Omega$
with inverse mass rescaling preserves the total conjugated
mass potential and leaves the curvature term unchanged.
Additional source or readout gain factors can absorb the
constant conformal field weights if absolute calibration is
released. Unspecified incoming signals can spoof a finite scalar
score. The declared retarded problem excludes them. Neither case
addresses a dynamical metric reconstruction theorem.

**Reviewer decision: B, effective theory.** The access-relative
two-site result and both weak scalar-channel distinctions survive
the independent checks. They follow from structural modification
and prescribed-metric Klein–Gordon propagation. The result
identifies what the extra interaction supplies, with no separate
GRUT dynamics or global geometry inference.
[Structural modification](https://doi.org/10.1006/jsvi.1999.2993) ·
[Conformal covariance](https://arxiv.org/abs/0806.0803) ·
[Prior feasibility review](FEASIBILITY.md)

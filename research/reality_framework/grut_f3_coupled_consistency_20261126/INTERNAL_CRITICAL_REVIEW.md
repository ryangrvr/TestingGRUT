# GRUT-F3 independent internal review

A separately tasked internal AI derived the reduction and constructed
an exact rational replay after freeze
0452465ea9a3ba9b29b2543395f66549331d295c. The replay imports no builder
code. External RAI review and experimental replication remain open.
[Frozen protocol](FROZEN_PROTOCOL.md) · [Parameters](PARAMETERS.json)

## Constraint closure and characteristic propagation

On a regular Legendre branch define the normal velocities $v^I$,
matter momenta $\Pi_I=\sqrt h\,P_{IJ}v^J$, and
${\cal H}_m=\Pi_Iv^I-\sqrt h P$. The derivative $P_{IJ}$ is
defined by contraction with the symmetric variation $\delta X^{IJ}$.
At fixed $\Pi,h$, varying the Legendre transform cancels derivatives
of $v^I$. Thus the momentum derivative of the matter Hamiltonian
smeared by $N$ is $Nv^I$, and its field derivative is

$$
-N\sqrt h P_{\phi^I}
-\partial_a(N\sqrt h\,P_{IJ}h^{ab}\partial_b\phi^J).
$$

The antisymmetric matter bracket consequently contracts one derivative
of $N$ or $M$ with $\Pi_I h^{ab}\partial_b\phi^I$. This is the
matter momentum constraint with the same metric structure function
as the pure Einstein bracket. The periodic spatial boundary removes
surface terms. Mixed normal brackets cancel because the matter
Hamiltonian is ultralocal in $h$ with no gravitational momentum
dependence and the Einstein variation with respect to the
gravitational momentum is proportional to the lapse. Spatial
deformations act by the Lie derivative. This closes the total
classical first-class algebra off shell on the regular branch.
[Teitelboim, hypersurface deformations](https://doi.org/10.1016/0003-4916(73)90096-1)

The conservation identity follows from
$T_{\mu\nu}=P g_{\mu\nu}+P_{IJ}\partial_\mu\phi^I\partial_\nu\phi^J$
and the scalar equations. Both identities hold for a broad class
of $P$. Neither fixes the scalar propagation symbol. On the chosen
rolling background the field-1 symbol is
$-\rho_X\,\omega_\mathrm{proper}^2+P_X|\mathbf k_\mathrm{proper}|^2$.
Both coefficients are positive for the branch. The tensor and the
constant spectator have the metric-null principal symbol. The nested
cones admit the same FLRW slices. Subluminality here is a difference
of propagation with a common foliation, consistent with covariant
stress conservation.

There are six spatial metric configurations and two scalar
configurations before removing the four first-class generators
and their gauge directions. This gives four physical local modes.
The harmonic gauge metric equations have the usual metric principal
part and the scalar kinetic Hessian enters the scalar principal part.
The reduced quadratic modes being quantized are the comoving curvature,
the separate spectator and the two transverse traceless polarizations.
On the rolling branch and with a finite nonzero-mode cutoff,
their positive kinetic terms provide ordinary positive-norm oscillator
quantization. The frequencies are time dependent. This construction
does not quantize the background, select a vacuum, or test the
operator constraint algebra beyond linear perturbation theory.
[Garriga and Mukhanov, perturbations](https://arxiv.org/abs/hep-th/9904176)

## Independently recovered background and finite-band gap

Energy conservation gives
$F=dx/d\ln a=-6x(1+2x)/(1+6x)$.
Matching $X_{\rm can}=XP_X$ and $U_{\rm can}=XP_X-P$ exactly
matches both the pressure and energy density. For $x>0$ the
canonical clock is monotone, so $U_{\rm can}$ can be defined along
it and the conservation law gives its background equation. The
production point $x=0$ represents algebraic limits of ratios.
The physical expanding rolling solution and comoving scalar
description require $x>0$.

The independent checker uses the logarithmic derivatives of the
active and canonical pumps. Their $x$-dependent terms are
$\tfrac12\ln[(1+6x)/(1+3x)]$ and
$\tfrac12\ln[(1+2x)/(1+3x)]$, respectively. With $D=d\ln(z/a)/dx$,
one has $W=1+FD$ and $W_x=F_xD+FD_x$.
It evaluates $A=W^2+(1-\epsilon_H)W+FW_x$.
On the whole interval $0<x\le1/16$, $|F|\le3/8$,
$|F_x|\le6$, and $|1-\epsilon_H|\le2$. The active log
derivative is positive and bounded by $3/2$, with derivative
magnitude at most $27/2$. The canonical log derivative is
negative and bounded in magnitude by $1/2$, with derivative
magnitude at most $5/2$. These bounds give
$$
|A_{\rm act}|\le1059/128,\qquad
|A_{\rm can}|\le1347/256,\qquad |A_{\rm metric}|\le1.
$$
They are within the preregistered caps without using an interpolation
of the four fixture points.

At $x=1/32$, the exact active squared speed is $17/19$.
The principal split from each luminal control is $2/19$.
The pump differences $A_{\rm act}-A_{\rm can}$ and
$A_{\rm act}-A_{\rm metric}$ are positive on this fixture,
with exact values $4498848/4561235$ and $112488048/159643225$.
Even treating their signs as unknown, the registered caps and
$H/k_\mathrm{physical}\le10^{-3}$ keep both normalized
squared-frequency splits above $1/10$. The curvature contribution
to the canonical difference is bounded by $9/500000$.
This is an instantaneous coefficient comparison in the declared
subhorizon band. Turning it into measured frequency data requires
an instrument and a background evolution estimate. The illustrative
precision reserve does not certify either.

~~~bash
python3 critical_replay_f3.py --output /tmp/grut-f3-critic.json
~~~

The checker evaluates all four frozen fractions, the pump and local
dispersion reserves, the rank-six coefficient system and the two
determinant-squared controls. It runs well below a second on the
review machine. A read-only comparison to the production archive
matched 107 common scalar entries exactly. The Hessian comparison
permutes the two independently chosen coordinate orders before
comparing entries. [Independent replay](critical_replay_f3.py)

## What imposing exact common propagation requires

For a metric-null covector let $u^I=k\cdot\nabla\phi^I$. In
coordinates $(X^{11},X^{22},X^{12})$, the symmetric kinetic
variation from polarization $\eta$ is proportional to
$$
S(u)\eta=(u^1\eta^1,\ u^2\eta^2,\ (u^2\eta^1+u^1\eta^2)/2).
$$
The nonmetric null symbol is minus $S(u)^{\mathsf T}H_X S(u)$.
With two independent gradients, choose a nonzero null vector
outside their span. The derivative of the map from the null
cone to $u$ has rank two there. If it were rank deficient,
a nonzero combination of the gradients would be proportional
to the normal of the cone. The image therefore contains an
open set in the $u$-plane. If the symbol annihilates every
polarization there, its coefficients, polynomial of degree
two in $u$, vanish everywhere. The directions $(1,0)$
and $(0,1)$ fix five independent symmetric Hessian entries.
The direction $(1,1)$ then fixes the cross entry between
$X^{11}$ and $X^{22}$. Exact coefficient elimination
confirms rank six.

This claim requires the whole polarization space in the
null kernel. Testing that the matrix has just one null
eigenvalue, or testing parallel gradients, is weaker. On each
connected open fiber of the realizable branch a zero full
kinetic Hessian gives the affine form $P=A_{IJ}(\phi)X^{IJ}
+ B(\phi)$. No condition on the remaining functions of the
fields follows from that algebraic principal-symbol test.

For the frozen determinant-squared term, the parallel-gradient
Hessian is nonzero but its projected null correction is zero.
Indeed the adjugate of the rank-one $X$ annihilates the gradient
direction, making the first determinant variation vanish. At
the independent-gradient fixture the null correction is the
diagonal matrix $\operatorname{diag}(-1/8192,-1/16384)$
in the $\Lambda=1$ convention. The latter test establishes
that the operator is not absent from the full principal theory.

## Matter dispersion and verdict

For the rotating complex scalar, the radial mass at fixed angular
rate $\mu$ is $m_r^2=2\lambda r_0^2$. The kinetic mixing at
quadratic order is $2\mu\,\delta r\,\partial_t(r_0\delta\theta)$.
The secular equation for positive frequencies is
$$
(\omega^2-k^2-m_r^2)(\omega^2-k^2)-4\mu^2\omega^2=0 .
$$
The gap is $m_{\rm gap}^2=m_r^2+4\mu^2$. Eliminating the radial
field at leading derivative order gives the frozen $P(X)$ with
$\Lambda^4=\lambda v^4$. Put $u=(k/m_{\rm gap})^2$ and
$\beta=4\mu^2/m_{\rm gap}^2=1-c_s^2$. The exact light-branch
correction to the low-frequency squared speed is
$$
\frac{\omega_-^2}{k^2}-c_s^2
=\frac{4\beta^2u}{(1+\sqrt{1+4\beta u})^2}\le\beta^2u.
$$
The benchmark bound is $1/9025$. This is a flat rotating-state
matching check, separate from the FLRW pump budget. A true
overlap band requires $H/m_{\rm gap}\le10^{-4}$ and controlled
derivative scales. The parameter hierarchy can be arranged,
but the combined curved-space radial corrections, loops and
experimental systematics have not been computed. The
microscopic canonical matter equations have metric-null fronts.
The reduced low-frequency cone and the canonical front address
different limits. [Joyce et al., relativistic superfluids](https://arxiv.org/abs/2204.03678)

**Reviewer verdict: C for the asserted common-cone implication.**
Identical homogeneous stress and metric evolution, common
first-class constraints and regular positive quadratic modes
have a counterexample with different effective cones. The
exact luminality classification shows what additional premise
forces an affine kinetic form in the stated class. The result
uses ordinary gravitational and matter effective theory.
No separate GRUT interaction, observational claim or full
quantum gravitational dynamics follow from this test.

# Principal-symbol obstruction and conditional restriction

## What the common source and constraints leave free

At a smooth background, let $q(k)=g^{\mu\nu}k_\mu k_\nu$
and $u^J(k)=\nabla^\mu\phi^J k_\mu$. For the scalar system
the principal matrix on a covector $k$ is
$$
 \mathcal P_{IK}(k)=P_{IK}q(k)-P_{IJ,KL}u^J(k)u^L(k),
 \qquad P_{IJ,KL}=\frac{\partial P_{IJ}}{\partial X^{KL}} .
$$
It follows by keeping the second derivative on the varied
scalar in $\nabla_\mu(P_{IJ}\nabla^\mu\phi^J)$.
Stress depends on at most first matter derivatives and the
scalar equations contain at most first metric derivatives.
Thus in hyperbolic metric gauge the second-order principal
blocks separate. The physical tensor symbol is metric-null.
An invertible local field-coordinate change transforms the
scalar principal matrix by a nonsingular congruence in its
principal part. Each physical cone is invariant under
coordinate or field relabeling.

On the explicit rolling solution the active mode obeys a
principal wave equation with proper-time kinetic coefficient
$D=P_{1,X}+2X P_{1,XX}=1+6x$ and gradient coefficient
$P_{1,X}=1+2x$. The spectator, tensor and matched canonical
active mode have the metric speed. The active speed is
$$
 c_{\rm a}^{\,2}=\frac{1+2x}{1+6x}.
$$
Both coefficients are positive and $c_{\rm a}<1$ for
$x>0$. The homogeneous time slices are spacelike for all
the symbols. More directly, on this branch the active symbol
is negative for every metric-timelike normal covector
because $P_{1,X}>0$ and $P_{1,XX}>0$. The effective signal
cone is nested inside the metric cone. These compatible
hyperbolic systems with positive kinetics establish the
obstruction to inferring a single metric characteristic
speed from common gravitational coupling, constraint
closure, conservation, or a common space of admissible
local clock slicings.
[Acoustic propagation](https://arxiv.org/abs/0708.0561)
· [Two-field characteristic examples](https://doi.org/10.1088/1475-7516/2008/04/017)

## Conditional classification of an exact common metric cone

Using a single metric as the characteristic standard for
all wave probes would impose a stronger requirement. Suppose **every
independent scalar polarization** lies in the kernel of
$\mathcal P(k)$ at each nonzero metric-null covector $k$.
At a fixed local derivative configuration with two linearly
independent scalar gradients, this implies
$P_{IJ,KL}=0$ there. Here two full polarizations, or a
two-dimensional kernel, is part of the hypothesis.
Coincidence of some roots of a determinant alone supplies
less information.

To prove the claim choose a null covector outside the
span of the two gradient covectors. Such covectors exist
in four-dimensional Lorentz space. Restrict the linear map
$k\mapsto (u^1,u^2)$ to the smooth nonzero null cone.
The derivative at the chosen $k$ has rank two on the
cone's tangent hyperplane. If it did not, a nonzero linear
combination of the scalar gradients would be normal to the
tangent hyperplane and proportional to $g^{-1}k$.
Lowering the index contradicts the choice of $k$.
Therefore the projected values $u$ contain an open subset
of $\mathbb R^2$. Since the full symbol vanishes for
these $k$, the quadratic polynomials
$P_{IJ,KL}u^J u^L$ vanish identically in $u$ for each
$I,K$.

Let $T_{ijkl}=P_{ij,kl}$ in field indices. It is symmetric
in the last pair. The vanishing quadratic says
$T_{ijkl}=-T_{ilkj}$. Chaining this relation with last-pair
symmetry gives
$$
 T_{ijkl}=-T_{iljk}=T_{iklj}=-T_{ijkl},
$$
so every component vanishes. The exact finite algebra
check evaluates the three independent entries of the
principal matrix on $(1,0),(0,1),(1,1)$ after extending the
polynomial. It has rank six on the six independent
symmetric-Hessian entries. These directions are test
vectors in the polynomial argument. Physical accessibility
of the projected directions used to prove openness is
established separately.

If the condition holds throughout a connected open,
realizable kinetic branch at fixed field values, the
vanishing Hessian integrates to
$$
 P(\phi,X)=h_{IJ}(\phi)X^{IJ}-U(\phi).
$$
On stable branches the arbitrary field-space metric
$h_{IJ}$ is positive definite. The converse follows from
$\mathcal P(k)=h\,q(k)$. The input of exact universal
metric principal propagation therefore excludes the
nonlinear first-derivative kinetic terms in this class
once imposed across that branch. It still leaves $h,U$
and the geometric action's coupling as free physical
inputs. Measuring one homogeneous trajectory, testing one
characteristic polarization, and imposing common
propagation only in a vacuum are weaker conditions.
The proof concerns first-derivative scalar actions with
metric dependence only through the displayed contraction.

## A homogeneous control with hidden nonlinear kinetics

Take the analytic comparison operator
$\Delta P=(\det X)^2/\Lambda^{12}$. On any homogeneous
two-field gradient both covectors are parallel and
$X=yy^{\mathsf T}/2$. The determinant vanishes and
$\operatorname{adj}(X)y=0$. Every first variation of $X$
at a homogeneous background, including from the metric,
is a sum of terms with at least one factor $y$.
Hence $\delta\det X=0$ as well. The added operator has
zero background value and zero quadratic action on these
backgrounds. It can be invisible even if both homogeneous
fields roll and the Hessian in the space of abstract
symmetric matrices is nonzero.

The registered gradients in units of $\Lambda^2$ are
$(1/4,0,0,0)$ for both fields in the parallel case, and
$(1/4,0,0,0),(0,1/4,0,0)$ in the independent case.
Take the metric-null covector $(1,0,1,0)$ in the
orthonormal chart. The respective contributions to the
principal matrix are
$$
 \Delta\mathcal P_{\parallel}=0,\qquad
 \Delta\mathcal P_{\rm independent}
  =\operatorname{diag}(-1/8192,-1/16384)
$$
in frozen dimensionless units. This is a local
principal-symbol control on declared derivative configurations.
The inhomogeneous fixture has no coupled cosmological
evolution assigned to it. [Exact fixture values](RESULTS.json)

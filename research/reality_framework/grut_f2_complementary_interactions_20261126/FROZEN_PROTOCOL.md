# GRUT-F2 frozen analytical protocol

26 November 2026. This packet continues from draft PR #24 at
f51fa0abbc046ea28b9e796fcad6e42d609ed86e. The frozen parameters refer to
the original F1 spring input by hash. F1's continuum comparison specifies
Minkowski versus a compact conformal bump outside the apparatus, without
a particular profile or finite observation window. The instances below
specialize that same class and retain the general F1 invariance.

## Question and data

Test whether the largest number of hidden sites jointly perturbed in a
setting is two in a worst-case reconstruction theorem for positive spring
networks measured at one unchanged port. In particular, test the sharper
hypothesis that finite diagonal loading patterns at one positive Laplace
parameter identify the graph even with dark modes, while the entire data
from every singleton pattern can fail. The intervention class consists of
static nonnegative independent diagonal loading at physically labelled
sites, with a reset between settings. The local actuator support is
calibrated without knowing spring connectivity; no hidden readout is
provided. Sequential within-trial changes and any other physical tests
are outside the lower bound.

For spacetime, test two independently calibrated Klein–Gordon probes
with different conformal transformation properties at the same apparatus
and finite times. Estimate an observable difference with explicit
propagation and error bounds and study a joint consistency relation.
Separate identifiability in a fixed model class from elimination of a
known rival and from global geometric reconstruction.

## Spring settings

The two F1 Hamiltonians have unit masses and common positive pinning,
with exactly their original positive edges. The response is the
root driving-point entry of $(V+zM+D)^{-1}$ with $z>0$ the square of
the Laplace variable, $D$ a static nonnegative pin increment, and $M$
the kinetic mass matrix. Record these controls.

1. With all interventions at the original port, retain the F1
   symmetry argument even for a nonlinear drive. Uniform mass doubling
   and a root pin shift are exact invariant numerical controls.
2. Use the independently tagged hidden site 3. At the frozen shifts
   evaluate a positive unit pin contrast and a unit inertial mass
   contrast at the original readout. Verify each distinguishes at
   $z=1$. At this shift the mass and pin matrices give the same scalar
   increment. Evaluate all three latent hidden tag choices at
   $z=1,2$ and keep collisions if any. The minimum extra contrast
   for the known pair is zero versus one when the calibrated site-3
   perturbation succeeds. No new site readout is used for this claim.
   Check the finite-time tail at cutoff 20 for the calibrated
   $z=1$ contrast using the positive spectra bound.
3. Recover the full stiffness of each F1 network from one root
   scalar per setting, at $z=1$, using no loads, two distinct
   positive singleton loads at each hidden tag, and one simultaneous
   positive equal load on each pair of hidden tags. Prove uniqueness
   for the class of connected positive equal-baseline-mass spring
   networks with a known positive pin, without assuming full
   Krylov rank. It is sufficient to prove a strictly monotone
   inverse on the physically admissible resolvent entries.
4. For the worst-case failure of singleton-only data, attach a
   root with equal positive springs to each vertex of the
   16-vertex rook graph or Shrikhande graph, both with unit
   baseline masses and pins. Verify the registered common
   adjacency polynomial and nonisomorphism. Check exact
   singleton port equivalence analytically for all positive
   Laplace shifts and strengths, and replay all singleton
   settings at two registered shifts. A double load on
   $(0,0)$ and $(0,2)$ must distinguish them.
5. Add an independent, observed unit-mass unit-pin oscillator
   $y$ with elastic energy $(y-b^{\mathsf T}x)^2/2$ in each network.
   A collective $b$ fixed by the F1 intertwiner is an invariant
   control. A channel supported on tag 3 must distinguish at the
   registered shifts. Evaluate the full two-port matrix.
6. A tagged hidden quartic $+\gamma x_3^4/4$ gives a Hamiltonian
   nonlinear control. Compute the leading power and coefficient
   of the term linear in $\gamma$ and cubic in the root kick,
   from a common corresponding harmonic preparation. Register
   distinct powers 21 and 13 on the cycle and pendant triangle.
   A quartic confined to the observed port, or a radial quartic
   on the hidden block, retains the intertwining.

## Spacetime settings

Use signature $(-,+,+,+)$ with $c=\hbar=1$, the no-boundary
globally hyperbolic Minkowski chart, and $\tilde g=\Omega^2\eta$,
$\Omega=1+\epsilon b$. The bump is centered at $(t,\mathbf x)
=(1,(1,0,0))$. Use the support box in time times the spatial ball
from PARAMETERS.json, with smooth nonnegative profiles and unit
interior plateau. Explicitly, convolve the even cutoff that equals
one for $|v|\le5/8$ and vanishes for $|v|\ge7/8$ with an even
normalized smooth mollifier supported on $[-1/8,1/8]$.
On the transition interval the cutoff is the descending quintic
$1-10u^3+15u^4-6u^5$ with $u=4(|v|-5/8)$. The mollifier is
proportional to $\exp[-1/(1-64v^2)]$ in the open interval and
zero elsewhere. Call the result $\zeta$ with support $[-1,1]$
and plateau $[-1/2,1/2]$. Set
$b=\zeta(8(t-1))\zeta(8|\mathbf x-(1,0,0)|)$.
The second derivative bound 240 and first derivative bound 8
imply the registered bound on $\Box b$.

The physical apparatus occupies the declared flat tube.
Use a smooth radial profile of integral one proportional to
$\zeta(8|\mathbf x|)$. The temporal envelope is one for
$|t|\le3/8$ and zero for $|t|\ge1/2$, smoothly and monotonically
joined by the cutoff $E(1-u)/(E(1-u)+E(u))$ with
$u=8(|t|-3/8)$ and $E(v)=e^{-1/v}$ for $v>0$ and zero otherwise.
Center the source envelope at zero and use the same envelope
at the detector times. Couple to the canonical scalar with the
same real strength and proper units in both geometries. Fix
$m$ and $\xi$ by independent calibration assumptions. In
particular flat-region calibration alone cannot determine
$\xi$. Use ordinary retarded conditions with zero prior source
response and no additional incoming boundary signal.

Compare the scaled real observable $Y=16\pi\int g_{\rm det}\phi\,d^4x$
from the retarded equation
$-(\Box_g-m^2-\xi R)\phi=f$. All support is inside the same
flat tube. The measurement interval ends before time 3.

* The massless conformal probe must be invariant exactly, for
  the whole restricted kernel.
* All fields are blind in the early gate by causal separation.
  A translated bump centered at time 4 is also blind in the
  whole finite measurement tube. Keep those negative controls.
* The massive conformal probe must give a signed late signal
  on the fixed weak and strong backgrounds. At the strong
  setting certify the stated strictly positive magnitude floor.
  The same positivity estimate tests strict monotonicity
  of the scalar score in the known profile family.
* The massless minimally coupled field must give a signed late
  signal for the same weak background, using a first Born
  formula and a uniform bound on the remainder. The
  strong-setting estimate for this field is a stress control
  and may remain inconclusive.
* Derive the rank-two compatibility equation for inferred
  effective potentials and a bound on the ratio of the two
  known-channel linearized score shifts. Calibrations,
  restricted data and any inverse step taken as an input must
  be declared. Compare with exact ambiguities under altered
  mass, gain or boundary assumptions.

## Gates and budget

One rational implementation and one independently authored check,
each within 30 seconds and 17 sites. Use analytical support results
and interval certificates in the continuum. No continuum grid or
parameter sweeps. A finite witness can verify the exact algebra,
with no tolerance or fitted rounded classification. Keep negative
and inconclusive controls. Verdict B fits successful conditional
reconstruction and same-apparatus discrimination explained by
ordinary structural modification and field propagation. Verdict A
requires a distinctive relationship that survives that comparison.
Residual ambiguities are reported at their own scope regardless of
the positive tests.

# What the joint physical response constrains

## Joint model and observable

The discrete map sends a fixed spring matrix in a tagged
site basis to the family of root responses indexed by
calibrated loads. Equality of measurements with supports
through size two forces equality of the matrix in the
declared class. The port-preserving transformations from F1
remain an equivalence if they also intertwine all supplied
load operators. Equal static singleton data can persist
even when the matrices are not related by a transformation
fixing all such operators simultaneously. The rook–Shrikhande
control illustrates that failure of joint isolation by
singleton data. The positive inverse uses coordinated inputs.
Unknown tag permutations can only rename the reconstructed
coordinates. No graph embedding into an ambient geometry
is assigned. [Spring theorem](SPRING_DERIVATIONS.md)

For the continuum, the common model class is the prescribed
four-dimensional Minkowski conformal class with positive smooth
$\Omega$ and $\Omega=1$ on a calibrated flat apparatus tube
$U_T$, and with separate Klein–Gordon species sharing this
metric. Fix the species masses, curvature couplings,
retarded boundary convention, and canonical source/readout
normalization. The output data for field $i$ are the restricted
operator
$$
\mathcal R_i(\Omega)= E_{U_T}^{*}\,
      (\partial_t^2-\Delta+m_i^2+V_i)_{\rm ret}^{-1}\,E_{U_T}.
$$
Here $E_{U_T}$ extends sources by zero. The tested records
are declared scalar pairings of each $\mathcal R_i$. The
coefficient of $\partial_t^2-\Delta$ and the dimension are
inputs. The three probe equations neither require a common
microscopic field nor change their background geometry.
An allowed diffeomorphism identity on $U_T$ acts on the
background and all couplings together, leaving the data.

| Property | Licensed conclusion |
| --- | --- |
| Causal structure | Shared causal support and the early blind interval follow the same principal cone and the access window. This gives consistency checks on the measured region. Here the global cone is assumed via the conformal class. |
| Conformal information | The massless conformal channel compares propagation in that class. Its equality cannot select $\Omega$. General causal-order or boundary reconstruction requires enough spacetime-resolved data and its own inverse hypotheses. |
| Metric scale | The massive channel can separate the two prescribed scales and identifies the amplitude in the known bump profile family. Local scales are calibrated on $U_T$. Arbitrary scale reconstruction requires additional inverse data. |
| Curvature | The minimal channel responds to $\Box\Omega/\Omega$, which is related to scalar curvature. The score is integrated and does not give a curvature map. |
| Dimension | Site count and the four-dimensional scalar equation are inputs. Neither joint test infers an ambient dimension. |
| Field parameters | Constant mass leaves a local flat tail signature. Flat local calibration leaves curvature coupling undetermined. Allowing unknown spacetime-dependent field potentials recreates a response ambiguity. |

## A joint compatibility equation

Put $A=\Omega^2-1$ and $B=\Box_\eta\Omega/\Omega$.
For each field the effective potential is
$V_i=a_i A+b_i B$, where $a_i=m_i^2$ and $b_i=1-6\xi_i$.
If two potential maps in the same chart have been obtained
by separate valid inverse measurements, and
$d=a_1b_2-a_2b_1\ne0$, they give
$$
A=(b_2 V_1-b_1 V_2)/d,\qquad
B=(a_1 V_2-a_2 V_1)/d .
$$
A common positive conformal scale implies $1+A>0$ and
the exact differential condition
$$
 B=\frac{\Box_\eta A}{2(1+A)}
   -\frac{\partial_\mu A\,\partial^\mu A}{4(1+A)^2}.
$$
With compatible support and boundary calibration, this is
also sufficient for those two zeroth-order coefficient
functions to arise from the claimed common $\Omega=\sqrt{1+A}$.
The equality follows from the chain rule and conformal
Klein–Gordon reduction. It tests the common-geometry hypothesis
when the individual potential maps and calibrations are
independent. Establishing those maps from a specific restricted
wave experiment remains a separate inverse problem.

For the frozen apertures there is also an operational
weak-response relation. Let $D_M,D_R$ be minus the derivative
at $\epsilon=0$ of the late scores in the massive conformal
and massless minimal channels respectively. The two are
positive for the chosen nonnegative nonzero $b$. With the
scaled flat waves $F_M=4\pi R_m f$ and
$W_M=4\pi R_m^{\rm adv}g_{\rm det}$,
$$
\frac{D_M}{D_R}
 =m^2\frac{\int_K b F_M W_M\,d^4x}{\int_K b|\mathbf x|^{-4}\,d^4x} .
$$
Both massive waves differ from $|\mathbf x|^{-1}$ by at
most $m^2/2$ on $K$. If $r_\pm=1\pm1/8$ and
$h=m^2r_+/2$, the ratio lies in
$$
[m^2 r_-^2(1-h)^2,\ m^2 r_+^2(1+h)^2]
 =[.0445461\ldots,.0847612\ldots].
$$
This assumes the same apparatus gains, and only the
nonnegative bump bounds, not its precise shape. It is a bound
on linearized responses on the prescribed family; the finite
signals have separately derived nonlinear error estimates.
A verified violation of the bound would reject the field,
geometry, aperture or calibration assumptions of this test.
[Continuum derivation](CONTINUUM_DERIVATIONS.md)

## Exact remaining equivalences

**Causal complement of the data.** Retarded operator
differences factor through the region in which the
coefficients differ. If no source-to-detector causal path
intersects it inside the window, both operators have the
same restricted data for every channel. The translated
bump is inequivalent curvature outside that reach. This
obstruction holds for the full joint operator on the
finite tube, even with calibrated parameters.

**Finite scalar records inside causal reach.** Two
independent scalar scores can distinguish a known rival
while still allowing exact coincidences for other profiles.
For the full smooth-profile class, take three linearly
independent nonnegative bumps supported inside the measured
causal diamond, each supported in the exterior of $U_T$.
Use independent small positive coefficients in an open
three-dimensional parameter box. Retarded scores depend
continuously on these coefficients by the same Volterra
bounds. The massless conformal channel is constant.
Two real scores cannot continuously inject an open
subset of $\mathbb R^3$ into $\mathbb R^2$, by invariance
of domain after adjoining a zero coordinate. These
profiles can be placed in the neighborhood $K$.
A conformal coordinate change fixing the tube pointwise
cannot rename profiles there: each point in that
neighborhood lies on at least two distinct straight
null lines meeting the open tube. A conformal diffeomorphism
of Minkowski fixes each such line as a set and fixes its
intersection with the other line. Thus the metrics
differ even modulo the apparatus-preserving diffeomorphisms.
The claim concerns finitely many measured scores.
It does not give indistinguishability for all of
$\mathcal R_i$ inside the reachable region.

**Unknown effective coefficients.** A constant mass on the flat
aperture can be fixed from the near-cone tail
$-m^2/(8\pi)$ compared with the local light-cone strength.
An independently fixed curvature coupling takes other
information because $R=0$ on the aperture. If a nominally
massless field may be conformally coupled, that field can
remain blind. Curvature amplitude and an unknown coupling
are degenerate at first order in a single weak contrast.
With known constant $m>0$, a curvature coupling cannot
cancel the scale term everywhere on the frozen profile:
the curvature contribution vanishes on the interior
plateau where $\Omega\ne1$. More generally, an exact global
coefficient cancellation by any nonzero fixed row
$(m^2,1-6\xi)$ would require
$(1-6\xi)\Box\Omega+m^2\Omega(\Omega^2-1)=0$.
Compact support in time and uniqueness of that semilinear
wave equation with flat past data force $\Omega=1$.
If $1-6\xi=0$ and $m>0$, the same result is algebraic.
This is an obstruction to hiding the bump at the full
operator-coefficient level. It implies no general
uniqueness for incomplete data with a free constant $\xi$.

In a larger class with unknown local scalar potentials,
any $\Omega$ can be compensated by subtracting
$V_i(\Omega)$ in the flat reduced operator, or equivalently
by supplying a covariant potential contribution
$-V_i(\Omega)/\Omega^2$ in the lower-order term on the curved
background. This is exact for every independently adjustable
channel and gives the same complete observed retarded map
even when the deformation lies within the causal diamond.
For the massive conformal probe the compensating effective
mass can be chosen positive, $m_i^2(x)=m_i^2/\Omega^2(x)$.
The free potential equivalence changes the calibrated model
class. [Conformal coefficient gauge in Lorentzian inverse problems](https://arxiv.org/abs/1607.08690)

**Gains and boundary conditions.** An unknown constant
gain can mimic a single massive score when both values have
the same nonzero sign. In the frozen late gate the flat
massive signal is a nonzero negative timelike tail, while
the curved massive score is still negative. Here
$m\sqrt{t^2-r^2}\le 3/4$ on the source-to-detector pairs,
where the alternating series for $J_1$ is positive. Ordinary local
calibration of the light-cone strength can fix the gain in
the product of source and detector. The late flat massless
signal is exactly zero by Huygens support, and a multiplicative
gain alone cannot fake the nonzero minimal-coupling contrast.
An unrelated source that contributes a common homogeneous
solution cancels when taking the response to the controlled
source. An unknown source-dependent boundary law can change
the inverse task. For example an interior domain with an
effective causal boundary impedance equal to the exterior's
retarded Dirichlet-to-Neumann map reproduces its interior
response by restriction and uniqueness. The frozen comparison
specifies the full-space retarded inverse, so this arbitrary
boundary law is outside its hypothesis.

The observed differences follow ordinary fields on prescribed
geometry. This hierarchy identifies where additional physical
information enters; it provides no dynamics for the metric,
coupling unification or quantum gravity.

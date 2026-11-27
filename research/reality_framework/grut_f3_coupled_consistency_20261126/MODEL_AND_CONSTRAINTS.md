# Dynamical metric, matter and canonical constraints

## Physical assumptions

| Item | Role in the calculation |
| --- | --- |
| Four dimensions and signature $(-,+,+,+)$ | Adopt a smooth Lorentzian spacetime and a spacelike foliation. Choose periodic spatial slices for the explicit example, so integrations by parts have no edge terms. |
| Metric and gravitational coupling | Vary one metric with the Einstein–Hilbert action and positive reduced Planck mass $M$. The example fixes the renormalized classical cosmological coefficient to zero. Universal sourcing by the variational stress and the field equations follow from this chosen action. |
| Matter representation | Adopt two real scalars with local first-derivative matter function $P(\phi^I,X^{IJ})$. Fix the regular Legendre branch. The sign and size of the kinetic correction are matter data. |
| Symmetry | Adopt diffeomorphism invariance. The active demonstration has a shift symmetry. The canonical comparator has a determined potential. |
| Clock and readouts | Along the representative rolling homogeneous solution use comoving proper time and its field value as a monotone local clock. Linear curvature on uniform active-field slices, transverse traceless tensors and the spectator are physical linear observables. |
| Quantum kinematics | Quantize the reduced linear modes with canonical commutation relations and a fixed cutoff regime. The clock choice, reference state and quantization are inputs. The linear response characteristics are independent of the Gaussian reference state. |
| Boundary and initial data | Set periodic boundary conditions, the expanding flat FLRW branch and constant spectator. Fix $a_*=1$ and the active kinetic density at the registered slice. Full state selection, cosmological history beyond the controlled branch and gravitational UV completion are open. |

Units set $c=\hbar=1$. The shared action is
$$
S[g,\phi]=\int \sqrt{-g}\,d^4x
 \left[\frac{M^2}{2}R+P(\phi^I,X^{IJ})\right],\qquad
X^{IJ}=-\frac12\,g^{\mu\nu}\partial_\mu\phi^I\partial_\nu\phi^J .
$$
The metric is dynamical and is varied in both comparison
models. Its history follows the field equations and initial data.

## Stress and matter symmetry

Define symmetric pair derivatives by
$\delta P=P_{,I}\delta\phi^I+P_{IJ}\delta X^{IJ}$ with full
indices summed. The Euler expressions and metric stress are
$$
 E_I=\nabla_\mu(P_{IJ}\nabla^\mu\phi^J)+P_{,I},\qquad
 T_{\mu\nu}=P_{IJ}\nabla_\mu\phi^I\nabla_\nu\phi^J+P g_{\mu\nu}.
$$
Direct differentiation and symmetry in the field labels give
$\nabla_\mu T^{\mu}{}_\nu=E_I\partial_\nu\phi^I$. This is
also the Noether identity from varying the matter action under
a diffeomorphism. On shell it gives covariant conservation,
consistent with $\nabla_\mu G^{\mu}{}_\nu=0$ and
$M^2G_{\mu\nu}=T_{\mu\nu}$. In the shift-symmetric active
example the separate field current is conserved as well.

## Normal and tangential constraints

Write the metric with lapse $N$, shift $N^i$, spatial
metric $h_{ij}$, and normalized scalar velocity
$y^I=(\dot\phi^I-N^i\partial_i\phi^I)/N$.
Then $X^{IJ}=(y^I y^J-h^{ij}\partial_i\phi^I\partial_j\phi^J)/2$.
Regularity means one can locally solve
$\pi_I=\sqrt h\,P_{IJ}y^J$ for $y$ at fixed $h$ and spatial
gradients. The matter constraint densities are
$$
\mathcal H_m=\pi_I y^I-\sqrt h P,\qquad
\mathcal D_{m,i}=\pi_I\partial_i\phi^I .
$$
With $\Pi^{ij}$ the conjugate momentum of $h_{ij}$, its trace
$\Pi=h_{ij}\Pi^{ij}$, and $R^{(3)}$ the spatial curvature,
the gravity densities are
$$
\mathcal H_g=\frac{2}{M^2\sqrt h}\left(\Pi^{ij}\Pi_{ij}-\frac12\Pi^2\right)
    -\frac{M^2\sqrt h}{2}\,R^{(3)},\qquad
\mathcal D_{g,i}=-2h_{ij}D_k\Pi^{jk}.
$$
The full Hamiltonian is $H[N]+D[\vec N]$ plus primary lapse
and shift momentum constraints. The smearings run over the
declared boundary class. [ADM](https://arxiv.org/abs/gr-qc/0405109)

For the matter-matter normal bracket, differentiating the
Legendre transform gives
$$
\frac{\partial\mathcal H_m}{\partial\pi_I}=y^I,\qquad
\frac{\partial\mathcal H_m}{\partial(\partial_i\phi^I)}
 =\sqrt h\,P_{IJ}h^{ij}\partial_j\phi^J .
$$
The terms involving $\delta y$ cancel at fixed momentum.
Only derivatives of the lapse smearings survive their
antisymmetrization. The combination
$y^I \sqrt h P_{IJ}=\pi_J$ then gives
$$
 \{H_m[N],H_m[L]\}
 =D_m[h^{ij}(N\partial_j L-L\partial_j N)].
$$
The matter density is algebraic in $h$ and independent of
$\Pi^{ij}$. The functional derivative of the gravity
normal Hamiltonian with respect to $\Pi^{ij}$ is local in
the gravity smearing. The mixed normal brackets therefore
cancel by antisymmetry in $N,L$. Adding the ADM brackets
and using the common spatial action on scalars, their
momentum densities, and the metric gives
$$
\begin{aligned}
\{H[N],H[L]\}&=D[h^{ij}(N\partial_j L-L\partial_j N)],\\
\{D[\vec u],H[N]\}&=H[\mathcal L_{\vec u}N],\\
\{D[\vec u],D[\vec v]\}&=D[[\vec u,\vec v]] .
\end{aligned}
$$
The same brackets hold for both actions, with the shared
metric structure function. The algebra is off shell on
their regular branches. Its evaluation at the matched
metric gives exactly the same structure functions.
The bracket depends on momentum transport and covariance.
The characteristic coefficients examined in this packet
remain independent data of the matter kinetic response.
[Hypersurface deformations](https://doi.org/10.1016/0003-4916(73)90096-1)

The reduced configuration count is $6+2-4=4$ local
physical modes after the four first-class constraints,
excluding the separate lapse and shift primary pairs.
They are two tensor and two scalar modes. Lapse and
shift enforce constraints and supply no extra kinetic
degrees. This is a classical gauge count.

## Exact dynamical background match

In the active model let $P_1(X)=X+X^2/\Lambda^4$ with
$X=X^{11}=\dot\phi^2/2$ and $P=P_1+X^{22}$.
The shift-symmetric spectator is constant on the
background. Put $x=X/\Lambda^4>0$. The exact densities and
equations are
$$
\begin{aligned}
 p/\Lambda^4&=x+x^2,&
 \rho/\Lambda^4&=x+3x^2,\\
 3M^2 H^2&=\rho,&
 \dot H&=-XP_{1,X}/M^2,\\
 a^3 P_{1,X}\dot\phi&=J,&
 \frac{dx}{d\ln a}&=-6x\frac{1+2x}{1+6x}.
\end{aligned}
$$
Choose $J$ fixed by $a=1, x=1/32$ on the expanding
$\dot\phi>0$ branch. The algebraic current relation and
Friedmann ODE define the homogeneous solution. They imply
$\dot\rho=-3H(\rho+p)$, and determine a monotone clock on
the rolling interval.

For the canonical comparison use scalar $\chi$ and potential
$U$, keeping the spectator and Einstein action. Define on
the same solution
$$
\dot\chi=\sqrt{2XP_{1,X}},\qquad
U(\chi(t))=XP_{1,X}-P_1=X^2/\Lambda^4 .
$$
Equivalently the local potential is given by the invertible
quadrature
$$
\frac{d\chi}{dx}=-\frac{M}{\sqrt6}
 \frac{1+6x}{x\sqrt{(1+2x)(1+3x)}},\qquad
U(\chi)=\Lambda^4 x(\chi)^2.
$$
This defines a covariant scalar potential on the interval
of interest. An additive choice of origin for $\chi$
renames its argument. Its kinetic energy is $XP_{1,X}$,
so its homogeneous $\rho,p$ match the whole original
history exactly. Conservation and $\dot\chi>0$ give its
canonical Klein–Gordon equation. Thus the same positively
expanding $a(t)$ solves both coupled systems, with the
same $T_{\mu\nu}(t)$, curvature, and initial metric data.
Their matter responses about this common history can still differ.
[Conventional equation-of-state degeneracy](https://www.sciencedirect.com/science/article/pii/S0927650507000904)

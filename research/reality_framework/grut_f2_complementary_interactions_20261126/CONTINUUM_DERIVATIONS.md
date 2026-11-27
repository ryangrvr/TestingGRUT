# Complementary scalar probes of the F1 spacetime comparison

## The invariant and the scale-sensitive potential

Use Minkowski $g=\eta$ and the F1 family
$\tilde g=\Omega^2\eta$, $\Omega=1+\epsilon b$.
F1 imposed a compact smooth bump with $\Omega=1$ in the
apparatus region. It left the profile free. The frozen F2
profile obeys exactly those restrictions. In particular it
has nonzero $\Box_\eta b$. Otherwise its vanishing data before
its compact support would force $b=0$ by wave uniqueness.
Thus $R_{\tilde g}=-6\Omega^{-3}\Box_\eta\Omega$ is nonzero
somewhere and this metric is inequivalent to Minkowski.
Both metrics have the same causal structure and agree on
all clock and length calibration in the aperture.
[F1](../grut_f1_operational_reconstruction_20261126/RESPONSE_AND_SPACETIME.md)

For each calibrated physical field let $u=\Omega\phi$.
The conformal operator identity is
$$
(\Box_{\tilde g}-R_{\tilde g}/6)(\Omega^{-1}u)
 =\Omega^{-3}\Box_\eta u .
$$
For the physical equation with source
$-(\Box_{\tilde g}-m^2-\xi R_{\tilde g})\phi=f$, the flat
reference equation becomes
$$
 \left[\partial_t^2-\Delta+m^2+V_{m,\xi}\right]u=\Omega^3 f,
 \qquad
 V_{m,\xi}=m^2(\Omega^2-1)+(1-6\xi)\frac{\Box_\eta\Omega}{\Omega}.
$$
On the source support $\Omega^3 f=f$ and on the detector
support $\phi=u$ and the integration volumes agree. Retarded
zero past data and finite propagation hold for the smooth
coefficients on both globally hyperbolic backgrounds.
The response in free quantum scalar theory has the same
kernel, since the field commutator is the causal propagator.
For $m=0$ and $\xi=1/6$, $V=0$ exactly, retaining the F1
negative control for the entire restricted kernel, at either
amplitude and with arbitrary source and detector in the flat
tube. [Conformal scalar law](https://doi.org/10.1007/s00220-009-0780-x)

## Finite observed interval

The apparatus lies in $(-1,3)\times B_{1/4}(\mathbf0)$ and
source and detector profiles have spatial support in
$B_{1/8}(\mathbf0)$. The bump region $K$ has time support
$[7/8,9/8]$ and spatial support in $B_{1/8}((1,0,0))$.
A source-to-bump or bump-to-detector spatial separation is
between $3/4$ and $5/4$. The source envelope is centered at 0
and the late detector envelope at 2. Their positive plateaus
contain the null-connected times of every point in $K$ for
every spatial point in the apparatus source support.
This is an aperture and pulse choice in the declared conformal
chart. It uses no calibration inside $K$.

Any first propagation from the bump to the apparatus support
can arrive no earlier than $7/8+3/4=13/8$. The early detector
closes at $3/2$. Its contrast therefore vanishes exactly for
all the calibrated fields, including their massive tails.
The translated deformation centered at time 4 is invisible
throughout the finite experiment by the same retarded support
argument. This second control holds for all retarded tests
with both endpoints in the full declared tube, not just
the two selected detector envelopes.

## Tail and repeated-scattering bounds

Write $R_m$ for the usual flat retarded inverse of
$\partial_t^2-\Delta+m^2$. Its kernel for $t>0$ is
$$
G_m(t,r)=\frac{\delta(t-r)}{4\pi r}
-\mathbf 1_{t\ge r}\frac{mJ_1(m\sqrt{t^2-r^2})}{4\pi\sqrt{t^2-r^2}}.
$$
One way to see the normalization is to use
$\rho=\sqrt{t^2-r^2}$ inside the light cone. The bounded
homogeneous radial solution is proportional to $J_1(m\rho)/\rho$
from the equation $h''+3h'/\rho+m^2 h=0$.
The wave operator on the tail step produces $4h(0)\delta(t^2-r^2)$
away from the tip, to cancel $m^2\delta(t^2-r^2)/(2\pi)$ from the
light-cone contribution. This gives $h(0)=-m^2/(8\pi)$.
The function $2J_1(v)/v$ is the average of $e^{iv y_1}$ on
the unit disk and has absolute value at most 1.
Thus the ordinary interior tail is bounded by $m^2/(8\pi)$.

For the forward flat wave $u_0=R_m f$ and the backward flat
wave $w_0=R_m^{\rm adv} g_{\rm det}$, the positive cone contribution
through the bump obeys
$$
4/5\ \le\ 4\pi R_0 f,\ 4\pi R_0^{\rm adv}g_{\rm det}\ \le\ 4/3 .
$$
Here the radial spatial profile has integral 1 and the retarded
or advanced arguments are on the unit plateaus. The absolute
massive change in either scaled wave is bounded by $m^2/2$
from its flat tail and the temporal support of length 1.
Inside $K$, for inputs supported there, both retarded and
advanced inverses have $L^\infty$ operator bound
$$
N_m=\frac{\Delta t^2}{2}+\frac{m^2\Delta t^4}{24}
=\frac1{32}+\frac{m^2}{6144},\qquad \Delta t=\frac14 .
$$
The spatial $L^1$ bound of the light-cone contribution at
time lag $v$ is $v$. The tail is bounded on the ball of
radius $v$ and adds at most $m^2 v^3/6$. Integrating these
bounds for $0\le v\le\Delta t$ gives $N_m$.

The change in the forward wave on $K$ solves
$u=u_0-R_m(Vu)$ there. If $\alpha=N_m\|V\|_\infty<1$, the
Neumann series bounds
$\|u-u_0\|_{\infty,K}\le \alpha\|u_0\|_{\infty,K}/(1-\alpha)$.
The same holds for the advanced wave. The measured late
signal difference for one field has the exact identity
$$
Y_\epsilon-Y_0
 =-16\pi\int_K w_0 V u\,d^4x,\qquad
Y=16\pi\int g_{\rm det}\phi\,d^4x .
$$
The overall factor is a fixed calibration in the chosen units.
The identity uses reciprocity of the flat retarded and advanced
inverses, valid for this real scalar equation.

## Massive conformal field

Set $m^2=1/16,\xi=1/6$. Then
$V=m^2(2\epsilon b+\epsilon^2 b^2)\ge0$.
For all $0\le\epsilon\le1/4$ the strong-setting bound is
$\alpha\le9219/8388608<.00110$. The free scaled wave is at
least $.76875$ through $K$. The perturbed scaled wave is
at least $.76724$ with the same conservative calculation.
Use the common floor $4\pi u,4\pi w>3/4$ on $K$.
The centered plateau $K_0$ has flat volume $\pi/24576$ and
$b=1$. Hence the exact signed late contrast satisfies
$$
 -(Y_\epsilon-Y_0)\ \ge\
 \frac{9}{16}\,m^2(2\epsilon+\epsilon^2)\,\frac1{24576}.
$$
At $\epsilon=1/4$ this is $27/33554432$ or approximately
$8.0466\times10^{-7}$. At $\epsilon=10^{-10}$ it exceeds
$2.8610\times10^{-16}$. These are guaranteed model signal
magnitudes, independent of numerical propagation.

More generally for $0\le\epsilon_1<\epsilon_2\le1/4$, insert
$V_{\epsilon_2}-V_{\epsilon_1}$ between the advanced wave
at $\epsilon_1$ and the retarded wave at $\epsilon_2$.
The same bounds give
$$
Y_{\epsilon_1}-Y_{\epsilon_2}
 \ge \frac{3}{1048576}(\epsilon_2-\epsilon_1)>0 .
$$
The single score thus identifies the amplitude in the supplied
one-profile family, with a quantitative slope bound. An unknown
bump shape is a different inverse problem.

## Massless minimal field

For $m=\xi=0$, $V_\epsilon=\epsilon \Box b/(1+\epsilon b)$,
where $\Box=\Box_\eta$. Its sign varies.
The two flat massless waves on $K$ are both exactly
$1/(4\pi|\mathbf x|)$. Their time envelopes are on the
plateau, and the integral of the radial normalized source
against the inverse-distance kernel equals the exterior
Coulomb potential by the spherical mean property.
Integration by parts using the compact support of $b$
then gives
$$
-\left.\frac{dY_\epsilon}{d\epsilon}\right|_0
 =\frac1\pi\int_K \frac{2b}{|\mathbf x|^4}\,d^4x
 \ge \frac{2(8/9)^4}{24576}.
$$
This uses $\Delta|\mathbf x|^{-2}=2|\mathbf x|^{-4}$ away
from the origin. Boundary terms vanish.

Let $B=32768$ bound $\|\Box b\|_{\infty,K}$. Indeed,
convolution preserves the normalized cutoff bounds
$\|\zeta'\|_\infty\le 8$ and $\|\zeta''\|_\infty\le240$.
Spatial radial differentiation contributes the additional
$2/r$ term only for $r\ge1/16$ from the bump center.
Together these give
$2\cdot 64\cdot240+2\cdot16\cdot8\cdot8=32768$.
For $\alpha_\epsilon=\epsilon B/32<1$, the error of replacing
the exact signal change by its first order is bounded by
$$
\frac{\operatorname{vol}(K)}{\pi}(4/3)^2\epsilon^2
 \left[B+\frac{B^2}{32(1-\alpha_\epsilon)}\right].
$$
The terms come from $V_\epsilon-\epsilon\Box b$ and from
the change of the forward wave respectively.
At $\epsilon=10^{-10}$ the remainder is less than 0.077
of the registered linear floor. The late negative
signal is therefore certified with magnitude greater than
$4.69\times10^{-15}$ in the same normalized units.
For the strong case the norm certificate is 256 and gives
no sign or accuracy guarantee. That stress case is inconclusive
under this bound.

These conditions distinguish nonzero responses within the same
measurement interval even though all probes have the same
principal causal cone. The long gate observes a tail sourced
through the bump. The weak-setting guarantees alone give no
evidence that a physical instrument can attain the very small
required error. [Rational certificates](RESULTS.json)

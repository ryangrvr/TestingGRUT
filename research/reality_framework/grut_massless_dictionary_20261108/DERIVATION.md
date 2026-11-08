# B1-01 — a finite-volume boundary and detector dictionary

This derives one explicit discretization of the declared experiment. It does not
claim that a finite spin chain uniquely determines physical continuum endpoints
or the subcell shape of a detector. Those choices must be specified. The new
coefficients below follow from geometry and Gauss law, not from minimizing any
massive or massless probability discrepancy. Prior M0/P1 results are already known.

## 1. Boundary, mass, length and time

Use hbar=c=g=1, even N, neutral sector, incoming link flux zero and theta=0.
The exact Gauss-law solution is
\[
L_j=\frac12\sum_{i=0}^j Z_i+\frac14[1+(-1)^j]
\equiv S_j+\frac14+\frac{(-1)^j}{4}.
\]
The smooth field is represented by S_j+1/4 after separating the explicit
alternating c-number. This does not assert that S_j has no further staggered or
finite-spacing fluctuations. The adopted continuum map has
Phi=phi/sqrt(pi)=E/g, Phi(0)=Phi(L)=1/4 in this boundary sector. Its basis and
limitations are given in [Okuda, equations (50), (62), (66)](https://arxiv.org/html/2210.00297v2).

The physical Hamiltonian is H/g=W_lat/(2sqrt(x)), x=(ag)^(-2). Set
m_lat/g=m_R/g-ag/8, with m_R=0 for the massless test. The bulk mass improvement is
the published convention of [Dempsey et al., equation (3), section IV.5](https://arxiv.org/html/2206.05308v3);
its periodic-translation derivation does not supply an exact open-boundary
matching theorem. No fitted mass, boundary term or energy-dependent clock is added.

Retain the previously declared ghost interval L=(N+1)a=5: sites
x_j=(j+1)a, links centered at (j+3/2)a. For the free hopping operator the exact
eigenfunctions sin[ell*pi*(j+1)/(N+1)] vanish at 0 and L, and
epsilon_ell=a^(-1)cos[ell*pi/(N+1)]. The half-filled gap is
2sin[pi/(2(N+1))]/a=pi/L-a^2*pi^3/(24L^3)+O(a^4/L^5).
This derives the free-theory ghost placement and velocity normalization.
Its use in the interacting finite chain is an inherited matching assumption;
we do not move the endpoints to improve the result. Time is the same physical
t in exp(-iHt) on both sides, not rescaled by a measured gap.

In the massless continuum, mu=g/sqrt(pi), and the quadratic boundary equation
(-partial_z^2+mu^2)bar(Phi)=0 gives
\[
\bar\Phi(z)=\frac{\cosh[\mu(z-L/2)]}{4\cosh(\mu L/2)}.
\]
The zero-boundary fluctuations have u_n=sqrt(2/L)sin(n*pi*z/L),
Omega_n=sqrt(mu^2+(n*pi/L)^2). The nonzero mean is retained. In the massive theory
this profile is only a basis displacement, not an imposed interacting mean.

## 2. Derive one finite-volume detector

The physical target is F=integral_0^L f(z)Phi(z) dz, with normalized Gaussian
f centered at L/2 of width sigma=0.8, c=2, W=exp(icF), K_s=(I+sW)/2.
f has dimensions inverse length; F and c are dimensionless. The instrument is
complete and its conditional state changes are part of the experiment.

Choose link cells I_j=[(j+1)a,(j+2)a], j=0,...,N-2. Define
\[
v_j=\int_{I_j} f(z)dz,\quad b=1-\sum_jv_j
=\int_0^a f+\int_{L-a}^L f.
\]
These coefficients are fixed error-function differences. The approximation
interprets S_j+1/4 as a cell average of Phi, and uses the fixed boundary value
on the two edge cells. It therefore specifies
\[
F_a^{\rm lat}=\sum_jv_j\left(L_j-\frac{(-1)^j}{4}\right)+\frac b4
=\frac14+\sum_jv_j S_j.
\]
The last equality is exact lattice algebra, not a state or mean fit. The
interior weights are not renormalized to one; b completes their partition of
unity. Compared with the old normalized midpoint-weighted raw detector, both
the smearing and its known c-number phase are changed by a declared rule.
This is a new approximation to the same target, not the original detector's
result relabeled. There is no free coefficient to tune.

Its matched continuum finite-volume observable is
\[
F_a^{\rm cell}=\int_0^L f_a(z)\Phi(z)dz+\frac b4,
\quad f_a(z)=\sum_j\frac{v_j}{a}\,\mathbf1_{I_j}(z).
\]
Compare lattice histories with both this cell reference and the ideal smooth
target. This explicitly separates the chosen spatial reconstruction from the
lattice-to-continuum identification. A cell-average interpretation is additional
finite-resolution dictionary data; it is not uniquely derivable from H alone.
Point evaluation of a continuum quantum field is not silently substituted for
a smeared observable.

## 3. Approximation order and its domain

For a classical C^1 profile y with the stated boundary values, integrating f
against cell averages and constant edge values has error bounded by
\[
\left|\int fy-\left(\sum_jv_j\overline y_{I_j}+b/4\right)\right|
\le a^2\|y'\|_\infty
\left[\frac{L-2a}{12}\|f'\|_\infty+\|f\|_\infty\right].
\]
For each cell use the covariance identity
Cov(f,y)=(2a^2)^(-1) integral integral [f(x)-f(z)][y(x)-y(z)] dx dz;
its bound is a^2||f'|| ||y'||/12. The edge bound follows from
|y(z)-y(0)|<=z||y'||, and similarly on the right. This proves O(a^2) for this
quadrature on the specified smooth class. Interpreting a midpoint value as a
cell average would add at most a^2||y''||/24 per unit total interior weight for
C^2 fields. Neither statement is a uniform operator bound on a quantum field.

For the massless quantum reference, compute the full Gaussian experiment for
f_a and f directly instead of applying a classical differentiability assumption
to vacuum fluctuations. The boundary substitution discards fluctuating edge
fields in F_a; their effect is included in the reference-to-target comparison.
No coefficient is selected from that comparison. The full lattice error can
still contain bulk dispersion, boundary matching and operator corrections;
we do not claim the whole dictionary has proven O(a^2) quantum accuracy.

## 4. Exact massless reference for the chosen cell smearing

For k_n=n*pi/L,
\[
f_{a,n}=\sqrt{2/L}\sum_j\frac{v_j}{a k_n}
[\cos(k_n(j+1)a)-\cos(k_n(j+2)a)],\quad
h_{a,n}=f_{a,n}/\sqrt{2\pi\Omega_n}.
\]
The mean is b/4 plus the integral of f_a bar(Phi), obtained by elementary sinh
differences. The ordered two-point kernel is
C_a(t)=sum_n h_{a,n}^2 exp(-i Omega_n t). Products of Weyl operators retain the
ordering phases from this complex kernel. First-readout probabilities obey
p_s=[1+s Re exp(ic mean-c^2 C_a(0)/2)]/2, while successive probabilities require
the full ordered bra/ket word sum. They are not products of first-readout laws.

The piecewise-constant f_a has total variation
B=|v_0/a|+sum_j |(v_{j+1}-v_j)/a|+|v_{N-2}/a|. Integration by parts gives
|f_{a,n}|<=sqrt(2/L)B/k_n, and therefore
\[
\sum_{n>M}h_{a,n}^2\le\frac{B^2L^2}{2\pi^4M^2}.
\]
Here B depends only on the registered smearing and cell geometry (for this
unimodal Gaussian, B<=2 max f). Use the inherited conservative ordered-history
tail bound 2^r r^2 c^2 times this variance tail for r readouts. This is a
massless Gaussian reference bound, evaluated in ordinary floating point; it
does not control a massive interacting parent. Mode doubling is reported too.

## 5. Disturbance and diagnostic decomposition

Use the same times (0), (0,1), (0,2), (0,1,2) and complete instruments in both
formulations. Report the first-readout TV separately from sequential TVs. For
each two-readout experiment, compare its final marginal with an isolated
readout in the stationary ground state. This measures disturbance from the
earlier unobserved readout. No branch is reset or conditioned away.

After the physical calculation, form diagnostic Gaussian summaries from the
lattice's own mean and ordered covariance C_lat(t), t=0,1,2. Positive spectral
weights of the ground-state covariance supply an auxiliary Gaussian oscillator
realization of these data; it is not the lattice state or a fitted dictionary.
Let G_lat use this complete two-point data; let G_static keep its mean and
variance but use normalized continuum-cell kernel C_a(t)/C_a(0). Then the signed
identity for each outcome distribution is
\[
p_{lat}-p_{target}=(p_{lat}-G_{lat})+(G_{lat}-G_{static})
+(G_{static}-p_{cell})+(p_{cell}-p_{target}).
\]
The four contrasts diagnose higher-than-two-point information, temporal
two-point response, static statistics, and spatial reconstruction respectively.
They are not additive causal error percentages; signs and cancellations matter.
These data-matched auxiliaries never replace the physical comparison and are
not transferred to the massive calculation. A remaining miss can be localized
in this diagnostic sense without assigning it uniquely to a microscopic
boundary operator or mass correction.

## 6. Conditional transfer

Only if every registered massless gate passes, use exactly the same v_j and b
rule at N20 and m_R/g=0.5. Keep its mass, Hamiltonian, preparation and time scale
unchanged. Compare with P1's fixed E14 smooth-target histories. No new bosonic
cutoff or Gaussian massive approximation is introduced. P1's interacting-tail
and overflow obligations remain unsatisfied even if this comparison improves.

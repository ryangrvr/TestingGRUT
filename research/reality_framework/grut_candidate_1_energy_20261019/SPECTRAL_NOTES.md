# Continuum determinant, cutoff bounds and numerical scope

The calculation uses the original profile, not a homogeneous approximation.
The original NPZ has SHA256
`0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3`.
The main subtraction/restoration definitions are in
[ENERGY_DERIVATION.md](ENERGY_DERIVATION.md). This note describes what is actually
integrated, and what the error statement does and does not mean.

## 1. Matrix continuum determinant

In the constant circular charged basis the free operator is diagonal, with
inverse lengths `kappa_+,kappa_-,sqrt(2+zeta^2)` as defined in the derivation.
The rotating local potential is the unitary transform of

\[
V=\begin{pmatrix}
S^2+.3F^2-1&0&2FS\\
0&S^2+.1F^2-1&0\\
2FS&0&3(S^2-1)+F^2
\end{pmatrix},\qquad V_D=.01(S^2-1).
\]

For regular columns `Y=rho^(ell+1) f`, solve

\[
f''+\frac{2(\ell+1)}\rho f'=(\operatorname{diag}\kappa_i^2+V)f.
\]

Use the regular power series through fourth order at the small origin radius,
then stabilize the stacked `(f,f')` by QR at intervals of one. Every triangular
determinant is retained. Match at R to the exact free decaying columns
`b_i=sqrt(rho) K_(ell+1/2)(kappa_i rho)`. The Wronskian matrix is proportional to
`f'+diag[(ell+1)/R-b_i'/b_i] f`; divide by the exact free regular normalization
`Gamma(ell+3/2)(2/kappa_i)^(ell+1/2)`. This produces the infinite-exterior Jost
determinant for a potential truncated at R, not a reflecting-box mode sum.

For its Born coefficients use the exact free half-line Green function. Write
`g_i(r)=r I_nu(kappa_i r) K_nu(kappa_i r)` and `e_i=b_i'/b_i`. The first
coefficient satisfies `B1'=sum_i V_ii g_i`. A Volterra matrix obeys

\[
H'_{ij}=(e_i+e_j)H_{ij}+g_i g_j V_{ji},\qquad
B'_2=-\sum_{ij}V_{ij}H_{ij}.
\]

These equations reproduce `Tr(G0 V)` and `-Tr(G0 V G0 V)/2`; they do not use
fitted weak-potential coefficients. The independent reporter block uses a
scalar Riccati equation, with its first and second perturbative equations
subtracted before integration. Its endpoint Wronskian correction is retained.

`spectral.py` preserves 2,362 executed point records, including vacuum,
weak-potential, radial-tolerance, endpoint and infrared controls. All sampled
Jost phases are close to zero. The code uses `slogdet` to avoid overflow; this
is not permission to discard a physical imaginary part. The root verifier
rejects a non-negligible phase. Positivity on this finite set of samples does
not prove the absence of every continuum instability.

## 2. Integration and independent checks

The first run integrates `ell=0..40` and `zeta=0..10` on intervals
`[0,.25,1,3,10]` with 6 and then 8 Gauss points per interval. It uses R=20,
origin radius .005, DOP853 tolerance `2e-7`. Symmetry logarithms are subtracted
and exactly restored, not deleted. The two finite-cutoff energy coefficients
are `-1.1000785657860743` and `-1.1000970384048012`.

The post-result audit is explicitly labeled as such in `audit_protocol.json`.
Its complete 8-point replay changes R to 25, the origin radius to .001 and
the tolerance to `1e-9`. It also calculates separate finite-difference block
determinants at `(ell,zeta)=(0,.3),(2,.7)` on grids of 250, 500 and 1000
interior points. Their Born coefficients come from direct tridiagonal inverse
matrices and traces, not the Volterra equations. The finest discrepancies in
the coupled Born-subtracted logarithms are approximately `1.59e-5` and
`7.74e-6`; each decreases by about four on halving the spacing.

The graph support is checked separately at R=20 and R=25 against the original
R=50 restoration. The final combination uses the R=25 graphs with the R=25
replay. Original results and the stricter replay are both retained; the latter
does not silently overwrite the first result.

## 3. Frequency tail: complex-resolvent-safe bound

Let `Vmax=sup ||V(r)||op`, including all four real species, and
`I2=int d^3r tr(V^2)`. Put `t=zeta^2+1-Omega^2`. The charged free resolvents
are `(h_l+t +/- 2i Omega zeta)^(-1)`. Heat-kernel positivity gives

\[
|G_{\pm,\ell}(r,r')|\le G_{t,\ell}(r,r')\le G_{0,\ell}(r,r').
\]

Thus one must not assume real Bessel arguments, but can still use positive
resolvent domination. In three dimensions

\[
\|G_0V\|\le V_{\max}/t=:r,\qquad
\|G_0V\|_{HS}^2\le I_2/(8\pi\sqrt t).
\]

For `r<1`, Schatten bounds and the logarithm series yield an order-three-and-up
remainder of at most `r ||G0V||HS^2/[3(1-r)]`. Therefore, at **all** ell,

\[
B_{\rm UV}=\frac{V_{\max}I_2}{48\pi^2}
\int_Z^\infty\frac{d\zeta}{(\zeta^2+1-\Omega^2-V_{\max})
\sqrt{\zeta^2+1-\Omega^2}}.
\]

With conservative `Vmax=2`, `I2=1800`, `Z=10`, this is
`0.03827594312628678` energy units `gv`. It is a bound on the subtracted
remainder, not an extra counterterm or a bound on higher loops.

## 4. Angular tail: two bounds, not a fitted extrapolation

The original compact-support estimate, checked in `CHALLENGE.md`, uses
`r_l=Vmax R^2/[(ell-1)(ell+2)]`. It is retained in `audit_results.json` and
gives conservative bounds .88915 (R20) and 1.98366 (R25) below Z=10.

The integrating root subsequently derives a sharper weighted bound. This is
an additional internal derivation, **not** a claim that the earlier challenger
audited it. Its protocol and outputs are `tail_refinement_protocol.json` and
`tail_refinement_results.json`.

Let `S=|V|^(1/2)`, `sigma=sign(V)` and
`T=S G0 S sigma`. Cyclicity gives `Tr T^n=Tr(G0 V)^n` wherever these traces
exist. The matrix-kernel norm is dominated by the scalar positive kernel
`sqrt(v(r)) g_l(r,s) sqrt(v(s))`, where `v(r)=||V(r)||op` and

\[
g_\ell(r,s)=\frac{r_<^{\ell+1}r_>^{-\ell}}{2\ell+1}.
\]

The radial centrifugal quadratic-form inequality
`h_l>=ell(ell+1)/r^2` gives

\[
\|T\|\le\frac{W_{\max}}{\ell(\ell+1)}=:r_\ell,
\qquad W_{\max}=\sup_r r^2\|V(r)\|_{op}.
\]

For the Hilbert-Schmidt norm, matrix Frobenius Cauchy-Schwarz bounds its kernel
by `g_l(r,s)^2 sqrt[tr(V(r)^2) tr(V(s)^2)]`. Using `2ab<=a^2+b^2` and the
exact integral over the squared massless kernel gives

\[
\|T\|_{HS}^2\le
\frac{2J_3}{(2\ell+1)(2\ell+3)(2\ell-1)},\qquad
J_3=\int_0^\infty dr\,r^3\operatorname{tr}V^2.
\]

The same logarithm bound now gives

\[
B_{\rm angular}\le\frac{ZJ_3}{3\pi}
\sum_{\ell>L}\frac{r_\ell}{(1-r_\ell)(2\ell+3)(2\ell-1)}.
\]

Use `Wmax=64`, `J3=600`, `L=40`, `Z=10`. The first omitted norm ratio is
`.0371661<1`; the bound is `0.05042781562786045`. The series is summed through
100000 and its remaining positive tail bounded by an `ell^-4` integral.
This estimate does not grow as R squared and also applies to a noncompact
profile **if** its Wmax and J3 are actually enclosed. It does not itself
enclose replacement of an exact soliton by the supplied numerical spline.

The angular estimate is used only below Z; the UV estimate covers everything
above Z. Their sum `.08870375875414723` has no double-counted region.

## 5. What is enclosed

The original profile is stored on cells of size .0025. `audit.py` converts each
local cubic for F and S to Bernstein coefficients. The convex hull, padded by
`1e-10`, encloses the interpolated values on each cell. Interval arithmetic
then encloses the symmetric local Hessian using a midpoint matrix norm plus a
radius row-sum bound. Upper rectangle integrals enclose its squared trace
moments. For R25 the computed bounds are approximately

| Input | Spline enclosure | Rounded input used |
|---|---:|---:|
| Vmax | 1.91284 | 2 |
| I2 | 1689.113 | 1800 |
| J3 | 550.727 | 600 |
| Wmax | 28.07922 | 64 |

Wmax is also checked on the entire supplied interval R50. These are padded
floating-point enclosures of an explicit cubic-spline model, **not** an
interval-arithmetic certificate for the exact boundary-value solution.

The original stationary solve, refinement and virial tests support that profile;
the domain/origin/tolerance replay supports the determinant numerics. Neither
is a rigorous nonlinear BVP error theorem. The final numerical range adds a
predeclared empirical `.02 gv` numerical allowance, a graph allowance and a
small classical/integer allowance to the analytical omitted-tail bounds, then
rounds the total **up**. It is not a statistical confidence interval, an
exact-profile theorem or an error bound on genuine two-loop, gravitational or
nonperturbative effects. Those qualifications accompany the binding result.

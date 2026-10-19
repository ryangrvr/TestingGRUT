# GRUT Candidate 1: fixed-charge one-loop energy and a free-carrier binding test

19 October 2026. New research calculation, continuing
[`1f370987f2fa96a52426ec95680a4f3ad7cbd373`](https://github.com/ryangrvr/TestingGRUT/commit/1f370987f2fa96a52426ec95680a4f3ad7cbd373).
The original Candidate 1 action, scalar couplings and MSbar matching are unchanged.

## Result

The actual spatially varying rotating object now has a numerical,
renormalized **one-scalar-loop fixed-charge energy calculation**, rather than
only a homogeneous effective-potential correction. At the predeclared nearest
integer charge `N=39944`, the calculation supports binding against breakup into
N infinitely separated elementary X carriers:

\[
\frac{K_1^{\rm ren}}{gv}\simeq 2.79528,\qquad
\frac{E_{\rm sol}^{0+1}(N)-N m_X^{0+1}}{gv}\simeq-3553.70573.
\]

The numerical allowance is **±0.11 in these energy units**, combining analytical
omitted-tail bounds for the supplied numerical profile with an explicitly
empirical numerical allowance. It is **not** a certified enclosure of the exact
continuum saddle, a statistical confidence interval, or a bound on genuine
higher-loop corrections. The fixed-charge semiclassical interpretation remains
conditional on stability of the reduced continuum saddle; finite-box spectra
and the current determinant samples support, but do not prove, that premise.

This is numerical one-loop support for one binding comparison. It does not
establish a quantum lifetime, all-channel quantum stability, empirical relevance,
or a universal theory of reality.

## 1. What was held fixed

Candidate 1 is Einstein gravity and supplied ordinary quantum matter plus a
complex charged scalar X, a real organization/binding scalar phi and a real
reporter D. Its new-sector tree potential is

\[
U=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2+g^2\phi^2|X|^2+
\lambda_X|X|^4+\frac12(m_{D0}^2+h^2\phi^2)D^2+
\frac{\lambda_D}{4}D^4.
\]

The original benchmark is `g=.1, lambda_phi=.01, lambda_X=.001, h=.01,
mD0=gv, lambda_D=.001, v/MP=1e-5`. Use the complete scalar potential and
counterterms required by the previous one-loop closure calculation. The
allowed additional scalar coefficients have their original zero finite
matching values at `mu0=gv`, not zero running at every scale. No finite
counterterm is fitted to improve binding.

This packet calculates flat-space **new-scalar loops**; SM and graviton loops
are outside scope. Spacetime, dimension, Einstein gravity, quantum probabilities
and the global U(1) structure are supplied. They have not been derived from
organization. A scalar vacuum-relative excitation energy does not set or
retune the total gravitating vacuum density.

The background is the original hashed `results_profile.npz`,
`sqrt(2) X/v=exp(-i Omega tau) F(rho)`, `phi/v=S(rho)`, `D=0`, where
`rho=gv r`, `tau=gv t`, and `Omega=.8000000019763136`. The prepared charged
matter is an initial resource; this is not formation from vacuum or a generic
autonomous-assembly demonstration.

## 2. Why this is a fixed-charge calculation

The rotating operator includes amplitude/phase velocity mixing, spatial
gradients, the organization field, and the reporter determinant. It is not a
local potential-curvature approximation. Its full definition and graph
restoration are in [ENERGY_DERIVATION.md](ENERGY_DERIVATION.md).

The unfixed-charge static amplitude curvature has a negative direction. We do
not replace that fact with a naive positive Gaussian. Routh reduction of the
parallel phase at fixed linear charge produces the positive rank-one term in
the **reduced** amplitude energy, while the original full gyroscopic operator
is unchanged. In an orthonormal kinetic basis,

\[
\det P_{\rm full}(\zeta)=\zeta^2\det P_{Q_1=0}(\zeta).
\]

The factor is retained when comparing with the complete original vacuum
determinant. An exact-rational canonical toy verifies the identity, the charge
constraint and equality between the logdet integral and the oscillator energy
difference. The latter agrees to `4.44e-16`. The toy also exhibits the negative
unfixed-charge curvature; it is a normalization check, not evidence of a
second physical soliton.

Phase and translation zero factors are integrable logarithms, subtracted only
for numerical integration and restored analytically with multiplicities one
and three. No vacuum oscillator or zero-mode finite energy is silently deleted.
The one-loop fixed-charge Legendre conversion is `E1(Q)=K1(Omega_class)`;
stationarity cancels the first-order frequency/profile-shift contribution to
the classical energy.

## 3. Actual continuum calculation

The numerical remainder is

\[
\frac1{2\pi}\int d\zeta\sum_\ell(2\ell+1)
\left[\log\frac{\det P_\ell}{\det P_{\ell,0}}-B_{1\ell}-B_{2\ell}\right].
\]

The determinant is computed from regular matrix radial solutions matched to
the exact decaying free solutions in the infinite exterior. First and second
Born terms use continuum free Green-function equations. MSbar-renormalized
one- and two-insertion graphs are then restored, including the actual
laboratory harmonics at Omega and 2 Omega. The reporter uses a separate scalar
Riccati determinant. [SPECTRAL_NOTES.md](SPECTRAL_NOTES.md) gives the method,
cutoff bounds and independent checks.

| Quantity, divided by gv | Value |
|---|---:|
| First restored graph, R25 | +2.83499381785 |
| Second restored graph, R25 | +1.05186785189 |
| Subtracted remainder, original R20, 6-point quadrature | -1.10007856579 |
| Subtracted remainder, original R20, 8-point quadrature | -1.10009703840 |
| Subtracted remainder, stricter R25 replay | -1.09158401905 |
| Selected finite-cutoff sum | +2.79527765069 |

The stricter replay increases R20 to R25, reduces the origin radius `.005` to
`.001`, and tightens DOP853 tolerance `2e-7` to `1e-9`. It re-executes all
`ell=0..40`, all 32 frequency nodes per ell, with no skipped points. Its change
`+.00851301935` is much larger than the initial 6/8 quadrature difference
`.00001847262`. Thus the initial quadrature agreement was **not** a complete
error estimate. Both original and stricter results are retained, and the
predeclared empirical `.02` numerical allowance exceeds the observed change.

Independent finite-difference determinants converge toward the continuum
answer at two nonsymmetry frequency points, with finest discrepancies in the
Born-subtracted logarithms of `1.59e-5` and `7.74e-6`. Their errors decrease
approximately fourfold when the grid spacing is halved. Their Born terms use
direct matrix inverses, not the same Volterra calculation.

All 2,362 original point records succeeded. Their determinant phases are below
`6.0e-14`; no non-negligible imaginary contribution was dropped. The stricter
replay also passes its phase gate. Fresh vacuum determinant controls at tighter
tolerance have errors below `8.7e-10`. These are finite numerical checks, not a
continuous-frequency stability theorem.

## 4. Omitted tails and numerical range

For the supplied spline, padded cellwise Bernstein ranges enclose the local
Hessian and its moments. Conservative inputs are `Vmax=2`,
`I2=int d^3r tr(V^2)<=1800`, `J3=int dr r^3 tr(V^2)<=600`, and
`Wmax=sup r^2||V(r)||<=64`. Actual spline enclosure estimates lie well inside
these inputs. They are not interval-certified errors on the exact nonlinear BVP.

The frequency tail uses complex-resolvent heat-kernel domination and a
Hilbert-Schmidt logarithm-series bound. The angular tail uses a weighted
Birman-Schwinger bound, not a fit to the last partial waves. Their regions do
not overlap: the angular estimate applies below frequency cutoff Z10, and the
frequency estimate covers every ell above Z10.

| Contribution to allowance, units gv | Size |
|---|---:|
| Analytical frequency tail, zeta>10 | 0.03827594 |
| Analytical angular tail, ell>40 below zeta=10 | 0.05042782 |
| Empirical numerical allowance | 0.02000000 |
| Restored-graph allowance | 0.00000100 |
| Classical/integer-state allowance | 0.00001000 |
| Sum, rounded upward for reporting | **0.11** |

The resulting numerical range for `K1/(gv)` is `[2.68528,2.90528]`. This range
is conditional on the stated spline/numerical approximations. Spatial endpoint,
origin and profile errors have numerical convergence evidence, not a rigorous
exact-saddle enclosure. Higher-loop, gravitational and nonperturbative errors
are not hidden in the `.11` allowance. Even the older, looser R25 angular bound
of about `1.984` would not change the free-carrier binding sign.

## 5. Integer charge and like-order binding

The original continuous classical branch point has
`Qref=39943.84624223011`, which is not an exact integer-charge eigenvalue.
`binding_protocol.json` selected `N=round(Qref)=39944` before the spectral result.
The strict energy conversion is

\[
\frac{E_{\rm sol}^{0+1}(N)}{gv}=
\frac{\bar E_0}{g^2}+\Omega_0(N-Q_{\rm ref})+\frac{K_1^{\rm ren}}{gv}.
\]

The integer correction is `+.12300621622`. A direct, additional fixed-N BVP
converges with maximum residual below `1e-10`; its charge is correct to the
quoted integration precision. It independently checks the energy expansion
to `2.13e-8` units gv. The sensitivity solve gives
`dQbar/dOmega=-4093.16324`, with an estimated quadratic integer correction
`-2.89e-8` units gv. The O(g^2) shift of K1 under this small change in Omega is
beyond strict one-loop order.

Use the previous packet's **strict first-order pole mass**
`mX^(0+1)/v=.10001232694040532`, not its curvature mass or a partial higher-order
square-root resummation.

| Same-order quantity, divided by gv | Value |
|---|---:|
| Classical energy at integer N, first-order interpolation | 36392.42286138 |
| Soliton energy including the calculated scalar loop | 36395.21813904 |
| N free-X pole-mass threshold | 39948.92387308 |
| Difference, soliton minus free carriers | **-3553.70573404** |
| Qualified numerical range for that difference | **[-3553.81574, -3553.59573]** |

The loop raises the classical energy by about `.007681%`. This is a different
calculation from the prior homogeneous bulk correction and must not be
identified with it. At this benchmark the energetic margin against complete
free-carrier dispersion remains large and negative at the calculated order.

## 6. Status after this calculation

| Claim | Current status |
|---|---|
| A specified scalar/Einstein EFT | Defined; supplied framework, not derived spacetime |
| Prepared classical charged organized object | Numerically supported in prior packets |
| Scalar one-loop counterterm closure and pole thresholds | Calculated in prior packets |
| Full rotating finite-charge one-loop energy | Numerically calculated here, with the qualifications above |
| Binding against N separated free X carriers | Numerically supported at one scalar loop |
| Residual soliton plus emitted X | Not calculated at one loop here |
| Multi-soliton quantum fission comparison | Not calculated at one loop here |
| Complete continuum stability theorem | Missing; finite numerical evidence remains |
| Quantum lifetime / tunneling rate | Not calculated |
| Autonomous formation from diffuse resources | Not established |
| Experimentally calibrated new sector / distinctive empirical detection | Not established |
| Universal theory of reality | Not achieved |

The achieved outcome is a **narrower working model with a demonstrated
classical regime and a numerical one-loop free-carrier binding calculation**,
not a validated or universal theory. The benchmark couplings are supplied
inputs, not a fit to observations. Known charge-stabilized scalar mechanisms
are not claimed as a new discovery merely by naming the sector GRUT.

## 7. Verification, provenance and next discriminating test

`VERIFICATION.json` records **38/38 checks passed**, including a fresh exact toy,
fresh continuum vacuum controls, an independent adaptive bubble-kernel
calculation, FD convergence, source/protocol identities, charge integrality
and energy arithmetic. The complete stricter spectral replay takes about
373 seconds; the original spectral pass about 291 seconds. Raw point records,
graphs, Fourier data, protocols and the original challenge file are included.

The constructor, calculator and challenger supplied complementary internal
work. The challenge file covers the formulas and initial controls it actually
reviewed; the integrating root reviewed the completed spectral output and
derived the sharper weighted tail bound. This is not an independent external
review, and worker agreement is not experimental validation.

The next calculation should compare **the same renormalized energy at several
integer charges**, especially N, N-1 and N/2, with uncertainties small enough
for the differences. That tests residual-object emission and the equal-split
channel, which one energy point cannot decide. See
[NEXT_CALCULATION.md](NEXT_CALCULATION.md). Do not change the action or fit a
finite counterterm to obtain a preferred sign. A failure is a result, not an
instruction to tune until it binds.

Reproduction: [README.md](README.md). Formulas:
[ENERGY_DERIVATION.md](ENERGY_DERIVATION.md). Numerical and bound details:
[SPECTRAL_NOTES.md](SPECTRAL_NOTES.md). Machine-readable integration:
[binding_results.json](binding_results.json).

# GRUT q-vacuum contribution: coupled growth, lensing and tensor response

21 October 2026. Codex researcher takeover and contribution.
Source: Claude's `8ec108c9807892a1c6e7afb81ce2b8a65864e4f1`.

**Result:** the chosen q action supplies a common response interface.
In its stated linear regime, growth, lensing and gravitational-wave
propagation are linked by the same K(q). A local representative has been
integrated with a dynamical background and controlled numerical error.
The small positive vacuum residual remains a supplied calibration in
the accelerating example. This work does not solve its selection.

This packet is additive. The mechanism is credited to Klinkhamer--Volovik;
its response limits are familiar scalar-tensor/f(R) physics. Our contribution
is their explicit application to the pinned GRUT construction, a worked
implementation and a common interface for checking several consequences.
Nothing is presented as a GRUT-specific new law or empirical confirmation.

## 1. The linked consequences

Let g_dyn denote the response of nonrelativistic matter, eta=Phi/Psi the
metric-potential slip, and Sigma the lensing response, normalized to the
same fixed G_ref. On the stable unscreened quasi-static branch,

\[
g_{\rm dyn}={M^2\over K}(1+s/3),\quad
\eta={1-s/3\over1+s/3},\quad
\Sigma={M^2\over K},\qquad
s={p^2\over p^2+K/(3K_R)},\quad p=k/a.
\]

For homogeneous geometric-optics propagation, after accounting for the
emission waveform,

\[
\left({d_L^{\rm GW}(z)\over d_L^{\rm EM}(z)}\right)^2
\simeq{\Sigma(z)\over\Sigma(0)}.
\]

The full derivation, potential-sign conventions, distinction between
curvature and trace masses, and regime assumptions are in
[DERIVATION.md](DERIVATION.md). The response functions cannot be changed
independently while retaining this action and regime. These relations
are shared with the matching f(R) theory and hence do not uniquely identify
GRUT. Screening at the source and detector can invalidate the use of
cosmological K in the displayed distance relation; no screened realization
is solved here.

Lensing response at fixed density perturbation differs from the final
lensing signal: the density perturbation itself grows differently.
Sigma near one therefore does not imply an unchanged lensing power
spectrum. No observed spectrum is calculated in this packet.

## 2. Worked local branch

Choose linear K(q) and quadratic U_mu(q), with positive kinetic coefficient
and mass. Eliminating q yields the exactly specified local representative

\[
S_g={M^2\over2}\int\sqrt{-g}\,[R+R^2/(6m^2)-2\Lambda].
\]

The positive constant Lambda is supplied in this example. In the parent
q description its value is part of the integration-constant/matching
choice. It is not derived by calculating matter growth.

The homogeneous field, Hubble rate, dust and radiation evolve by the
coupled field equations. The GR expansion history is not imposed on this
solution. Approximate subhorizon growth is integrated alongside it.
All inputs, initial conditions, units and calibration status are in
[INPUTS.json](INPUTS.json).
Dust and radiation are supplied matter descriptions. Their existence,
particle content, primordial abundances and quantum probabilities are
not derived by the vacuum-response calculation.

The supplied reference densities are Omega_m=0.3, Omega_r=0.00009,
Omega_Lambda=0.69991; a runs from 0.1 to 1. The following masses are
**numerical hierarchy controls**, far below the assumed unscreened
laboratory floor when H_ref is identified with the observed Hubble scale.
They are not observationally viable unscreened benchmark fits.

| m/H_ref | maximum absolute H/H_GR minus 1 | final D/D_GR, k/H_ref=30 | final D/D_GR, k/H_ref=300 |
|---:|---:|---:|---:|
| 300 | 0.00375793 | 1.06351010 | 1.33322463 |
| 600 | 0.00093612 | 1.02060708 | 1.27248378 |
| 1200 | 0.00023382 | 1.00560876 | 1.17896939 |

The matched initial growth values are D=a and dD/dln(a)=D at a=0.1.
This is a stated finite-interval preparation, not a calculated primordial
growing mode. At m/H_ref=300 and k/H_ref=300, the endpoint gives
g_dyn=1.16662317, eta=0.71428998 and Sigma=0.99996520.
The background can remain close to GR while a short-wavelength growth
response is appreciable. This is a known f(R) mechanism, reproduced here
for the explicit local branch.

The Lambda=0 control at m/H_ref=1200 decelerates at the endpoint:
q_dec=0.48591823, with q_dec=-1-Hdot/H². Thus this low-curvature,
near-equilibrium, zero-residual example does not produce current cosmic
acceleration. It does not exclude every initial condition or every
high-curvature history of q-theory.

## 3. Numerical and approximation checks

- The full Friedmann constraint stays below 4.2e-15 relative to the
  stated sum-of-absolute-terms denominator in the sampled background runs.
  The exact equations propagate its residual as Cdot=-4HC.
- Tightening rtol from 1e-9 to 1e-11 and atol from 1e-12 to 1e-14 changes
  sampled H by at most 2.1e-11 relatively and D by less than 4.2e-14.
  These are numerical comparisons for the specified equations.
- The exact constant-field de Sitter solution is retained to floating-point
  precision. The separate Einstein--de Sitter growth controls reproduce
  the analytic exponents 1 and (sqrt(33)-1)/4 to about 2e-13 in endpoint D.
- A direct complex-wave integration on a separately prescribed smooth
  background agrees with the tensor endpoint law within 8.2e-7 in the
  three controls; errors are not monotonically decreasing with frequency.
  That background is not asserted to solve the q cosmology.

There are eight numerical dispositions in RESULTS.json. They are not
eight independent empirical tests. The tight integration error does not
certify the quasi-static approximation: the smallest sampled p/H is
17.21 for k/H_ref=30 and 172.10 for k/H_ref=300. Maximum H²/(p²+nu²)
is 0.001686 and 0.00003343 respectively in the m=300 run. These are
scale estimates; no total relativistic perturbation-error bound was proved.

Full scalar perturbations, independent scalar isocurvature, nonlinear
screening, an observational likelihood, transitions and radiative
stability are not computed here. See SOURCES.md for exact read levels.

## 4. A conditional scale calculation

If a uniform unscreened scalar has m>=4 meV in both the laboratory and
the cosmological background, its Compton wavelength is at most
4.9332e-5 m. The additional force relative to the lensing normalization
at k=1 Mpc^-1 is then, in the flat/heavy approximation,

| Redshift | g_dyn/Sigma minus 1 |
|---:|---:|
| 0 | <=8.5198e-55 |
| 1000 | <=8.5369e-49 |

The 4 meV number is an assumed control relayed from Claude's ledger;
it is not an independently fitted bound here. These estimates concern
the scalar-mediated enhancement, not an arbitrary free drift in K or
every screened model. In the simple quadratic branch the adiabatic drift
is also suppressed by H²/m². A light cosmological mode can be considered
only with a specified environment-dependent model or another explicit
change; it is not obtained by silently dropping local constraints.

## 5. What this establishes for construction

We have a calculable common response, a solved local background example,
and several consequences at fixed inputs. This is useful structure for a
working hypothesis. It leaves the main dependency unchanged:

**What covariant dynamics or independently justified state rule fixes the
positive total renormalized residual, including its sign and magnitude?**

If a later version supplies such a mechanism, apply this interface to its
actual solution and determine which relations survive. A changing mu or
additional operator may alter the equations and perturbations; its source
stress and energy transfer have to be calculated as part of that version.
If the residual is calibrated instead, retain that label and test a
different consequence. Candidate 1 stays paused under its existing status.

## Reproduction and development record

Run from this directory:

```bash
python calculate.py
```

Dependencies are numpy and scipy; exact environment versions and source
hashes are in [RESULTS.json](RESULTS.json). No packages were installed.
The initial execution completed the integrations but failed to serialize
a NumPy boolean in the numerical-check dictionary. The serializer now
casts those dispositions to Python bool; the complete calculation was
rerun successfully. No physical parameter changed for that repair.
There is no preregistration claim, held-out comparison, outside review,
or independent laboratory measurement in this contribution.

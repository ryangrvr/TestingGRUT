# Challenger audit: rotating fluctuations and vacuum pole threshold

19 October 2026. This audit retains Candidate 1's action, original matching, and classical profile. It covers the actual fluctuation calculations as they arrive and an independent one-loop vacuum X pole-sign check. It does not calculate a finite-charge quantum energy or permit promoting a classical mode calculation into such an energy.

## Immediate interpretation gates

- A finite-box list of real frequencies is not a continuum completeness theorem. The box discretizes continuum scattering states and lifts zero modes; refinement and volume behavior must be reported.
- A negative eigenvalue of a static local Hessian is not itself a growing mode of the rotating coupled pencil.
- The full gyroscopic pencil uses the unmodified grand-canonical spatial operators. Adding the fixed-Q energy rank-one correction to Lu while retaining the original phase dynamics would double count the constraint.
- A pole mass calculated in the vacuum updates a physical free-particle threshold. It does not update the soliton energy or establish quantum binding until the finite-charge renormalized correction is also calculated.
- No sum of unweighted rotating-frame box frequencies is a renormalized fixed-Q one-loop energy.

## Direct pole sign and normalization check

At the tree vacuum the new scalar masses are `MX^2=c v^2=0.01 v^2`, `Mphi^2=2 lambda_phi v^2=0.02 v^2`, and the cubic tensor element is `g3=lambda_phi_x_x=2 c v=0.02 v`. Define the finite MSbar functions

\[
A(t)=t[\ln(t/\mu^2)-1],\qquad
B_0(s;a,b)=-\int_0^1 dz\,\ln\frac{za+(1-z)b-z(1-z)s-i0}{\mu^2}.
\]

Use the Euclidean inverse two-point convention

\[
\Gamma_E^{(2)}(p_E)=p_E^2+M_X^2+\Pi_E(p_E).
\]

The scalar one-loop effective action is `Gamma1=(1/2) Tr log H`. Expanding to second order in an external x background gives two off-diagonal `g3*x` insertions. The quadratic contribution is `-(g3^2/2) x I_E x`, so the inverse-propagator bubble is `Pi_E=-g3^2 I_E`. The two internal fields are different (phi and x), and the ordered tensor sum cancels the usual one-half: there is no additional one-half in this mixed bubble.

After continuation `p_E^2=-s`, the renormalized inverse propagator contains `-g3^2 B0(s)/(16 pi^2)`. Consequently

\[
\boxed{m_{X,\mathrm{pole}}^2=m_{X,\mathrm{curv}}^2
-\frac{(2cv)^2}{16\pi^2}\,[B_0(M_X^2;M_\phi^2,M_X^2)-B_0(0;M_\phi^2,M_X^2)]}
\]

through strict one-loop order, with the previously calculated tadpole-induced shift already included in the curvature. A positive sign in front of this standard B0 difference is wrong in this convention. Changing the definition of the self-energy or B0 can move signs between formulas, but must reproduce this inverse-propagator/pole relation.

There is an independent zero-momentum normalization check. The direct X curvature correction before shifting v is

\[
\frac{cA(M_\phi^2)+4\lambda_XA(M_X^2)+\alpha A(M_D^2)
-(2cv)^2B_0(0;M_\phi^2,M_X^2)}{16\pi^2}.
\]

At matching alpha=0. Here `B0(0)=1-2 ln 2`, and the expression reproduces the earlier CW X curvature correction `+3.36558452599e-7 v^2`. The `4 lambda_X` tadpole coefficient includes both real X modes. The distinct tadpole-displacement contribution is `2 c v delta_v`; it must not be omitted or counted again when starting from the strict shifted curvature.

At these exact tree masses the momentum-dependent difference is analytic:

\[
B_0(M_X^2)-B_0(0)=1+\ln2-\pi/2=0.122350853765\ldots
\]

and therefore

\[
(m_{X,\mathrm{pole}}^2-m_{X,\mathrm{curv}}^2)/v^2
=-3.09918333078\times10^{-7}.
\]

Using the previous strict shifted curvature gives `mX,pole^2/v^2=0.0100024653880811`. The square root of this truncated squared result is `0.100012326180732`; the strict first-order expansion of the mass itself is `0.100012326940405`. Their difference is induced higher-order content and must not be described as a computed two-loop correction.

The mixed scalar cut begins at `(Mphi+MX)^2`, about `0.05828 v^2`, above the on-shell `s=0.01 v^2`. Thus this one-loop bubble has no absorptive part at the evaluated pole. That statement concerns this vacuum scalar-loop channel, not all nonperturbative or gravitational channels.

Primary confirmation was inspected in S. P. Martin, *Two-loop scalar self-energies in a general renormalizable theory at leading order in gauge couplings*, [hep-ph/0312092](https://arxiv.org/pdf/hep-ph/0312092), equations (1.4), (3.1), and (3.3): the pole equation is `m^2-s+Pi(s)=0`, the scalar cubic contraction has its ordered-index one-half, and the scalar bubble function is `B_SS=-B`. No two-loop result is imported here.

## Rotating operator and fixed charge

For

\[
X=e^{-i\omega t}[f+u+iw]/\sqrt2,\qquad\phi=\varphi+\eta,
\]

the real equations contain `+2 omega wdot` in the u equation and `-2 omega udot` in the w equation. With `exp(-i nu t)`, the off-diagonal pencil entries are `-2 i omega nu` and `+2 i omega nu`. The spatial operators are those derived directly from the action, with `-omega^2` in both charged components and no such term in the neutral eta or D operator.

The exact charge expanded to quadratic order in this parametrization is

\[
Q=\int[\omega((f+u)^2+w^2)-(f+u)\dot w+w\dot u]d^3x.
\]

Its linear part is `delta Q=integral(2 omega f u-f wdot)`. All nonzero-frequency eigenmodes of a correctly solved full pencil obey delta Q=0 by charge conservation. The zero/generalized-zero sector must be treated separately; numerical projection of an arbitrary amplitude vector is not a substitute for the phase-momentum constraint.

At fixed linear charge, eliminating the component of wdot parallel to f produces the radial energy-Hessian correction

\[
L_u\longrightarrow L_u+\frac{4\omega^2}{I}|f\rangle\langle f|,
\qquad I=\int f^2d^3x.
\]

This belongs to a reduced quadratic-energy argument. It is not an extra potential term to insert into the unmodified full u/w gyroscopic pencil. A correct calculation must state which representation it evolves or diagonalizes.

## Zero modes, thresholds, and all-angular-momentum arguments

The global phase mode is `w=f` at zero frequency. The three translations have `u=partial_i f`, `eta=partial_i varphi`, `w=0`; they lie in ell=1. A spherical solution has no additional independent rotational zero modes. A time translation is proportional to its internal phase shift, not a new independent collective coordinate here. Their finite-box lifts must converge toward zero instead of being silently discarded or interpreted as a physical instability.

In the asymptotic vacuum the charged branches are `sqrt(k^2+MX,vac^2) +/- omega`; the full pencil also contains their negative-frequency partners. A Dirichlet radial box discretizes k. Its first positive scattering frequency need not equal the continuum edge, and a frequency below the finite-box level is not automatically a bound mode below the true continuum threshold. Partial-wave thresholds and localized eigenfunctions must be compared with volume dependence, including possible resonances. The neutral phi and D branches do not receive the charged rotating-frame shift.

An all-ell conserved-energy argument can be stronger than a finite spectral scan, but its premises must actually be established. The w operator has a nodeless zero mode f in ell=0; the eta/u translation pair supplies an ell=1 zero mode. Increasing ell adds the positive centrifugal form `[ell(ell+1)-2]/r^2` to the latter block for ell>=2. To infer positivity from the translation mode, one must establish that it is the lowest coupled ell=1 eigenfunction, for example using monotonic profiles and the appropriate cooperative-operator theorem after a sign change in one component. A computed near-zero translation residual alone does not prove that no lower eigenfunction exists.

The ell=0 amplitude/organization block still needs its fixed-Q correction and positivity check. D positivity follows directly at the matching point from its positive mass and nonnegative couplings. Conserved nonnegative quadratic energy can exclude finite-energy exponentially growing linear modes only after the physical constraint and symmetry directions are handled. It does not establish nonlinear orbital stability, arbitrary-data survival, or a quantum lifetime. Finite discretization evidence for the premises must not be relabeled a rigorous continuum theorem.

## Actual source and output audit

The audit now includes `DERIVATION.md`, `modes.py`, both numerical amendments, all 40 spectra in `mode_results.json`, the pole and full-tensor neutral-pole sources/results, and the angular-comparison source/protocol/results. `MODE_NOTES.md` and `POLE_THRESHOLDS.md` were checked against these sources and raw results. No new physical evolution or spectral run was performed by this audit.

The rotating generator has the correct gyroscopic signs and factors. The code uses the unmodified Lu in its dynamical matrix, and puts the fixed-charge rank-one term only in a separate quadratic-energy diagnostic. Its discrete charge-left-null vector, phase vector and differentiated-background Jordan partner agree with the analytic equations above. The constructor's canonical and velocity representations are consistent: their different-looking quadratic Hamiltonians are related by `r=qdot+Omega Jq`, rather than by a change of action. The dimensionless action carries the original `g^-2` normalization; a future quantum normalization must retain it.

The constructor's phase and coupled translation ground-state identities also pass the algebraic check. In the coupled ell=1 block, the change of sign in u makes the off-diagonal coupling negative. Under the stated strict monotonicity and regularity hypotheses, the positive translation pair then gives the displayed nonnegative form. This is a sufficient conditional theorem. A numerical translation residual, numerical monotonicity, and finite-grid radial inertia do not certify its continuum hypotheses. Likewise, the radial slope criterion additionally assumes an invertible amplitude operator with exactly one negative eigenvalue; a negative charge slope alone does not prove stability for an arbitrary branch.

The consistent backgrounds are solved at **fixed Omega**, not at exactly identical charge on every grid. Their Qbar values are 398.596525626, 399.062504484, 399.361905993 and 399.419277406 for N=120,180,400,800 at L=30, converging toward the original continuum reference 399.438462422. The fixed-charge statement is about perturbations of each such discrete saddle. It must not be described as exact matching of the prepared charge across all discretizations. The matched-spacing L30/L40 control removes this principal mesh confounding in the box comparison.

### Finite spectra and zero modes

All 40 declared spectra completed, with no skipped grid or recorded ARPACK nonconvergence. The dense matrices are fully diagonalized at the declared coarse grids. The fine grids use the declared bounded sparse shifts; they are not full fine-grid diagonalizations. The actual ell=5 amendment includes a complete 1080-root dense calculation at N180,L30 and three sparse controls.

No nonsymmetry internal frequency was found in the declared window `1e-5 < nu < .199998`, ell=0,...,5. No nonsymmetry growing root was found among the returned roots. The largest computed absolute nonsymmetry real part is `1.9224e-14`; the largest generator residual is `1.5656e-10`. The former is not an upper bound on a missed continuum growth rate. Nonzero-mode charge overlaps reach `7.2675e-13`, supporting the implemented charge identity without substituting for spectral completeness.

The initially positive phase root `s=+3.4022e-7` at N180 was preserved. Tightening the stationary solve changes it into the imaginary near-zero pair while reducing the exact phase/charge residuals to about `7e-15` and `4e-15`. Its phase overlap is effectively one. This supports identifying that pair as a numerical splitting of the defective symmetry zero. It does not authorize deleting a future converged positive real root merely because it is small or has a substantial symmetry overlap.

The translation frequency decreases from `.00704572` to `.00106142` over N120 to N800 at L30. At matched h=.05 it is `.00141706334` in both L30 and L40 boxes. This supports a discretization lift of the translation zero, rather than a volume-independent physical oscillator. The original roots remain in the raw data.

At matched h=.05, the first nonsymmetry positive frequencies are:

| ell | L30 | L40 |
|---:|---:|---:|
| 0 | .20661186569 | .20354848812 |
| 1 | .21171086850 | .20643457717 |
| 2 | .21666015239 | .21017881045 |
| 3 | .22674909993 | .21514328231 |
| 4 | .23652700375 | .22070862612 |
| 5 | .24749970926 | .22698922830 |

Their substantial volume dependence and approach toward the charged edge `.199999998024` preclude calling these six numbers internal excitation gaps. The numerical vacuum check independently reconstructs the correct discrete charged branches and the eta/D branches. The static charged tail scale `sqrt(1-Omega^2)≈.6` is not the rotating continuum edge.

The separate finite-grid fixed-charge radial amplitude minimum is positive, `.366966510879` at h=.05,L40, while the uncorrected grand-canonical amplitude eigenvalue is negative, `-.188833915084`. This is the expected distinction between the constrained energy diagnostic and an unconstrained static direction. It does not turn the negative static direction into a calculated instability. The constrained minimum lies near the static essential threshold `.36`; it is not a dynamical mode frequency or a new internal gap. At the same grid, ell=2 amplitude and phase minima `.285521978057` and `.380710532650` are positive. The centrifugal ordering for higher ell is exact for those discrete operators; promoting its numerical base positivity to a continuum assertion remains conditional.

### Angular comparison

The angular diagnostic correctly subtracts the asymptotic vacuum local matrix. Its perturbation matrix in (u,w,eta) is

\[
\Delta V=\begin{pmatrix}
S^2+.3F^2-1&0&2FS\\
0&S^2+.1F^2-1&0\\
2FS&0&3(S^2-1)+F^2
\end{pmatrix}.
\]

For real `|nu|<1-Omega`, the constant vacuum pencil is positive. If a genuine entire-half-line bound `Delta V >= -B/rho^2` holds, adding the centrifugal and radial kinetic forms therefore excludes a kernel whenever `ell(ell+1)>B`. This proves the stated conditional real-subgap exclusion. It is distinct from the conserved-energy argument excluding exponentially growing modes.

The source obtains the correct least eigenvalue of Delta V. Its optimized numerical maximum is `B=28.0100296087` near rho=4.2459486063 on the saved domain 0<=rho<=50. The comparison value 31 is a chosen numerical margin, not a certified upper bound on the exact half-line supremum. The code explicitly sets `half_line_bound_certified=false`. Small sampled tail values and a vacuum-valued endpoint do not establish the exact uncomputed tail bound or control every spline/background error. Thus ell>=6 exclusion remains conditional. Filling ell=5 closes the declared finite numerical angular workload, but does not remove that condition.

Even a certified version would exclude real modes only in the **common lowest gap**, `0<|nu|<1-Omega`. It would not exclude higher-channel bound states, embedded states or resonances. D has an independently positive operator and lies above this common gap; D states below its own, larger channel threshold are a different question. None of these arguments constitutes a complete continuous-spectrum construction.

### Vacuum pole calculations

The actual X-specific implementation agrees with the independent inverse-propagator derivation and analytic B0 difference above. Its tadpole displacement is counted once, its momentum-dependent mixed bubble has the correct negative sign, and its strict shell is evaluated at the tree squared mass. Its ten declared consistency checks pass. Solving the truncated inverse propagator implicitly or taking a square root without first-order expansion introduces higher-order content; the report correctly labels both separately.

The neutral-species route contracts the full real-field quartic and cubic tensors. Its quartic term is `0.5 sum K_iikk A(mk^2)`, and its bubble term is `-0.5 sum g_ikl^2 B0(s;mk^2,ml^2)/(16 pi^2)`, including ordered pairs. This correctly distinguishes identical-particle and mixed-particle loops. The separate vacuum displacement is `sum g_iia delta_v_a`. All four zero-momentum curvatures reproduce the preceding effective-potential result, x/y degeneracy holds, and the X pole agrees with the independent mixed-bubble calculation. Its six checks pass.

| Field | Strict one-loop pole squared / v^2 |
|---|---:|
| phi | .0200095687556280 |
| x and y | .0100024653880811 |
| D | .0100998306010239 |

All nonzero cubic channels are below their pair cuts at the evaluated external tree shells. This is a real, positive vacuum scalar-pole result within the pure-new-scalar one-loop approximation. It is not a soliton decay-width or lifetime calculation. The X pole-squared scale spread `1.94e-9 v^2` is a residual higher-order diagnostic, not a rigorous error estimate.

Replacing the thresholds in the uncorrected classical rotating operator by these vacuum pole masses would mix orders. A consistent loop-level soliton comparison also needs the shifted background/frequency, derivative and nonlocal fluctuation terms, renormalization, charge treatment and finite-charge energy. The earlier homogeneous CW potential has a negative local Hessian in part of the wall; integrating its clipped real part would not supply those missing contributions. The vacuum poles do not repair that omission.

## Bounded rejection tests and final verdict

The following would require rejecting a stronger claim:

- A nonsymmetry root with positive real part that persists under stationary-tolerance, radial-spacing and box controls would reject the corresponding no-growth claim. A classification threshold cannot override those controls.
- A controlled negative direction of the fixed-charge energy form would reject the claimed local energy-minimum argument. On its own it would still require analysis of the full gyroscopic dynamics before being called exponential growth.
- A localized, volume-stable real mode in the lowest gap would reject an asserted absence of internal modes and would have to enter any later determinant. The present empty finite-window lists do not rule out an arbitrarily shallow state.
- A frequency sum without continuum normalization, collective coordinates, ultraviolet subtraction/counterterms and the fixed-charge conversion would fail as a quantum energy, whatever its numerical sign.
- Comparing a corrected particle pole with an uncorrected classical soliton energy would fail as a quantum binding test.

**Verdict:** no sign, normalization or tensor-multiplicity blocker was found in the audited calculation. The outputs support the narrower statement that the declared tree scalar finite-box searches found no nonsymmetry growth or internal mode in their stated window, together with positive discrete constrained-energy diagnostics and conditional analytic routes to broader exclusion. The vacuum scalar pole correction is separately calculated and cross-checked. The audit rejects a complete continuum spectrum, a certified all-ell stability theorem, nonlinear stability, a renormalized finite-Q quantum energy, quantum binding or a lifetime as conclusions of this packet. Gravity, Standard Model loops and all ultraviolet completion claims remain outside this calculation. No finite-Q quantum energy has been calculated by this audit.

## Integrating report review

The numerical values and operator signs in `REPORT.md` were checked against the actual outputs. The report correctly distinguishes fixed-charge perturbations of discrete fixed-Omega backgrounds, classical versus loop-corrected thresholds, numerical growth searches versus continuum certification, and vacuum pole masses versus the missing finite-charge quantum correction. Two wording corrections were returned to the integrating worker: call this execution bounded parts of the preceding two-stage plan, since the requested continuum boundary conditions remain unimplemented, and say that the derivation proves the conditional implication rather than proving that its continuum hypotheses hold for this profile. These corrections change no calculation. The main claim is support for the declared bounded tests, with the stable-quantum-configuration milestone still open.

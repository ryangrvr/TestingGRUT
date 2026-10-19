# Fixed-matching scalar one-loop potential

The scalar Coleman–Weinberg diagnostic is complete and frozen. All physical matching values remain unchanged. The source is `effective_potential.py`; reproduce from any working directory with its full path, or from this directory with:

```sh
python effective_potential.py
```

The numerical protocol was written before execution to `potential_protocol.json`. Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 were used. The final run completed in 1.58 seconds; no installation, paid compute, or background calculation remains.

## Conventions and complete scalar Hessian

Canonical real fields are q=(phi,x,y,D), with X=(x+i y)/sqrt(2), and the evaluated constant backgrounds are (phi,F,0,0). In units v=1 the fixed potential is

    U0 = C + mphi2*phi²/2 + lambda_phi*phi⁴/4
           + mX2*(x²+y²)/2 + lambda_X*(x²+y²)²/4
           + c*phi²*(x²+y²)/2
           + mD2*D²/2 + lambda_D*D⁴/4
           + d*phi²*D²/2 + alpha*D²*(x²+y²)/2.

The last term is alpha*D²|X|². This alpha normalization agrees with the independent RG protocol. The matching values are lambda_phi=.01, lambda_X=.001, c=.01, d=.0001, lambda_D=.001, alpha=0, mphi2=-.01, mX2=0, mD2=.01, C=.0025, at mu=.1. A matched alpha=0 does not assert that its beta function vanishes.

At y=D=0, the full four-real-field tree Hessian consists of the block

    H_phi,phi = mphi2 + 3*lambda_phi*phi² + c*F²,
    H_x,x     = mX2 + c*phi² + 3*lambda_X*F²,
    H_phi,x   = 2*c*phi*F,

and the two remaining eigenvalues

    H_y,y = mX2 + c*phi² + lambda_X*F²,
    H_D,D = mD2 + d*phi² + alpha*F².

The code constructs the complete quartic tensor in all four real fields before differentiating. It therefore retains the off-diagonal derivatives needed for D and y curvature corrections, even though y=D=0 at the evaluated background. Both real X modes are counted once; no extra multiplicity is added.

For tree Hessian eigenvalues t_a=m_a²,

    V1 = sum_a t_a² [log(t_a/mu²) - 3/2] / (64*pi²).

MSbar subtraction and the specified finite constant are used, without fitting any new finite counterterm. For t<0 the code retains the complex principal logarithm log(t)=log|t|+i*pi, flags the negative Hessian mode, and reports the imaginary part. No negative mass is clipped. The continuous t² log(t) limit is used only for the potential at t=0; real analytic curvature routines require strictly positive tree masses.

## Vacuum, tadpole and constant-field curvatures

At the tree point (1,0,0,0), U0=0. The one-loop tadpole is

    dV1/dphi = -2.43874796154e-6 v³.

The strict order-hbar shift is delta_phi/v=+1.21937398077e-4. The shifted local vacuum is therefore phi/v=1.000121937398 at strict one-loop order. Directly minimizing the truncated function U0+V1 gives phi/v=1.000121881804, with residual gradient 2.24e-17 in the chosen units. The difference contains induced higher-order terms from minimizing a truncated potential; it is not an independently calculated two-loop correction.

| Field | Tree curvature m²/v² | Strict one-loop curvature m²/v² | Minimum-of-truncated-potential curvature m²/v² |
|---|---:|---:|---:|
| phi | 0.020000000000 | 0.020012778428 | 0.020012783213 |
| x and y, each | 0.010000000000 | 0.010002775306 | 0.010002775710 |
| D | 0.010100000000 | 0.010099830632 | 0.010099830626 |

These are zero-momentum constant-field curvature quantities, not pole masses. Their relative shifts are about +0.0639%, +0.0278%, and -0.00168%, respectively. All four local vacuum curvatures remain positive. This establishes a local vacuum check, not a global search of every vacuum.

The strict expression includes both the direct CW Hessian at the tree point and the tree-Hessian change induced by delta_phi. The direct X curvature correction alone is +3.36558452599e-7 v²; omitting the tadpole-induced vacuum shift would miss most of its total strict correction.

## Homogeneous depleted-phi bulk diagnostic

The requested local diagnostic is

    h(phi,F) = 2 [U0(phi,F)+V1(phi,F)-Veff(vacuum)] / F².

At tree level its depleted-phi minimum is phi=0, F=10^(1/4)=1.778279410039, with h0=.003162277660168 v². The four sorted tree Hessian masses squared at this point are

    (.003162277660168, .009486832980505, .010000000000000, .021622776601684) v².

They are all positive. The numerical region was frozen as phi in [-.15,.15], F in [1.4,2.2]. The tree Hessian is positive throughout this rectangle: conservative block bounds give A>=.0096, B>=.00588, |offdiag|<=.0066, hence determinant>=1.2888e-5; the transverse-X and D bounds are .00196 and .01. Corrected curvatures and the diagnostic Hessian were also checked on the frozen 31-by-101 grid. Their sampled minimum eigenvalues were .00195916 v² and .00220455, respectively. These grid checks are not a formal interval proof for the loop-corrected entire rectangle.

The strict order-hbar results are

    delta h = +1.17872091832e-7 v²,
    delta h/h0 = +3.72744282759e-5  (about +0.00372744%),
    delta F/v = +1.86702353203e-4,
    h_strict = .003162395532260 v².

Minimizing the truncated U0+V1 expression instead gives phi=0, F/v=1.778466020398 and h=.003162395556574 v². A one-dimensional minimizer and three independent two-dimensional starts agree on the minimum value to about 4.4e-19 v² and on F to about 2e-10 v. No point in this local rectangle has an imaginary CW term.

The square-root diagnostic divided by the tree X mass is .562351807730. Dividing instead by the truncated vacuum X curvature gives .562273777704. Both are homogeneous potential proxies. A curvature mass is not a pole mass, and this comparison does not establish a corrected physical free-carrier threshold or a quantum Q-ball binding energy.

## Negative-Hessian regions and the classical wall

Of 7381 points in the frozen background rectangle phi in [0,1.2], F in [0,2.4], 4342 have at least one negative tree-Hessian eigenvalue. This fraction describes the chosen numerical grid; it is not a probability or a claim that 58.8% of physical states are unstable. The largest principal-branch |Im V1| on that grid is 4.97359197162e-7 v⁴. All signs and complex values are retained in `potential_background_scan.csv`.

As a limitation check, constant-background Hessians were evaluated at field values along the original classical Q-ball profile, on rho=0,.02,...,20. Eighty-four sampled points have a negative mode, spanning sampled rho=2.64 through 4.30. The minimum squared mass is -.00357857280730 v² and the maximum local |Im V1| is 6.36927306327e-8 v⁴. Endpoints are limited by the .02 sampling interval. The raw flags are in `potential_profile_flags.csv`.

Those flags concern homogeneous static backgrounds with the same field values. They do not include the spatial-gradient operator, the rotating phase, the charge constraint, or the coupled fluctuation problem of the finite soliton. They neither contradict the previous constrained classical stability checks nor calculate its decay rate. The local complex potential was not integrated over the profile and was not treated as a real quantum soliton-energy correction.

## Vacuum energy and gravity calibration

With the fixed C=.0025 v⁴, the strict one-loop scalar vacuum energy is

    Vvac_strict = -1.22652632365e-6 v⁴.

The minimum of the truncated U0+V1 potential has Vvac=-1.22667495222e-6 v⁴; the difference -1.48628564202e-10 v⁴ is induced higher-order content, not a two-loop calculation. Because the tree vacuum energy is zero, a relative loop correction to that zero is not defined. The required strict constant shift is about 0.0490611% of the specified tree constant C.

To set the scalar contribution to a zero gravitating vacuum at strict one-loop order would require an additional matching/calibration shift +1.22652632365e-6 v⁴ in the total renormalized vacuum constant. An equivalent shift could reside in the independent gravitational vacuum parameter. No such shift was fitted or applied here. For a nonzero target total vacuum energy, the required shift changes by that target value and by any other sectors' vacuum contributions.

At the illustrative v/Mpl_reduced=1e-5, the raw strict scalar shift is -1.22652632365e-26 Mpl_reduced⁴. If no other gravity vacuum constant cancels it, it induces Lambda/Mpl_reduced²=-1.22652632365e-26 and R/Mpl_reduced²=-4.90610529461e-26. Its AdS curvature radius would satisfy L_AdS*v=1.56394757366e8. These are conditional gravity consequences of the uncancelled scalar shift, not a matched prediction of the observed cosmological constant.

Subtracting Veff(vacuum) in h(phi,F) defines a relative matter-energy reference only. It does not remove the gravitating vacuum energy, supply its physical calibration, or solve the cosmological-constant problem.

## Numerical verification, correction record and scope

Analytic spectral/Frechet derivatives include the degenerate x/y modes. Central-difference checks at three four-field backgrounds used steps 1e-3,3e-4,1e-4,3e-5,1e-5. At step1e-4 the maximum CW gradient and Hessian-from-values discrepancies were 1.05e-13 and 1.07e-13. Hessians differentiated from analytic gradients reached discrepancies of order 2.5e-15 at step1e-5. Values-based second derivatives show the expected cancellation loss at the smallest steps. The independent tree-potential second-difference check verifies quartic normalization. All raw differences are in `potential_results.json`.

Review caught a bookkeeping error in the first output: `tree_background_masses_squared` reused the mass array at the optimized truncated bulk point. The source now keeps `mbtree` separate and labels the optimized-background masses explicitly. This changed that metadata field only; no tadpole, potential, minimum, curvature correction or gravity result changed. The original output/log are preserved as `potential_initial_results.json` and `potential_initial_execution.jsonl`, and `potential_correction_record.json` records the correction. The final rerun passed the same frozen protocol.

The calculation demonstrates small scalar constant-field one-loop shifts at the fixed MSbar matching point and makes the static potential's limitations explicit. It omits derivative terms for the rotating phase, field-dependent kinetic corrections, nonlocal fluctuation effects, the finite-charge determinant, nonperturbative effects, and quantum lifetimes. It is not a full soliton one-loop calculation. No new operator was fitted, and the original candidate and classical-stability files were not modified.

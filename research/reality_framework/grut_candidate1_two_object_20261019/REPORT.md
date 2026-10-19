# Candidate 1: prepared objects interact under one action

19 October 2026 · Author numerical research result · Independent scientific review pending

**Result:** the same Candidate 1 scalar equations support several different, resolved finite-time two-object encounters. At the tested initial quarter-cycle phase, the two incoming charge concentrations interact and emerge as two separated outgoing concentrations with unequal charge sharing. Equal phase produces one central axial density component at the endpoint. Opposite and nearby phases produce two separated components. These qualitative distinctions survive the three tested resolutions.

This is progress from a prepared isolated object to a prepared interacting regime. It is not generic formation, an asymptotic scattering theorem, a long-lived many-object phase, a memory apparatus, a quantum calculation, or evidence that these fields exist in nature. In particular, “central component” does not certify a final stable fused soliton, and two axial regions do not alone certify elementary particles.

## Branch, frozen inputs, and execution evidence

- Research branch: [codex/grut-candidate1-two-object-20261019](https://github.com/ryangrvr/TestingGRUT/tree/codex/grut-candidate1-two-object-20261019).
- Parent: [96a7e7e](https://github.com/ryangrvr/TestingGRUT/commit/96a7e7e6a49f49dc30e1c92f0855522af41d9641).
- Kernel, profile, driver and 27-case protocol frozen before the full evolutions: [52f5528](https://github.com/ryangrvr/TestingGRUT/commit/52f552802759392a2353af068f60efea2c5c95b2).
- Completed suite and raw diagnostic outputs: [879dbd5](https://github.com/ryangrvr/TestingGRUT/commit/879dbd59cc307d03f8f2cf23c63bc7ad2ba00aac).
- [Completed GitHub run, including execution and verification logs](https://github.com/ryangrvr/TestingGRUT/actions/runs/39336072781).
- [Protocol](protocol.json), [machine-readable summary](results/summary.json), [execution amendment](EXECUTION_AMENDMENT.md), [coarse V8 records](results_v8/), and [supplementary controls](results/additional_controls.json).

The ordinary shell workspace was unavailable. Eight coarse evolutions were first executed in the available V8 runtime. The unchanged kernel, input profile, driver and physical protocol were then run in a branch-scoped standard GitHub-hosted job using Node v22.23.3. All 27 prescribed cases completed. The 33 numerical consistency conditions passed. All 27 output files were subsequently retrieved and the summary checks recomputed.

The eight coarse cases were compared between runtimes across 49,424 numerical entries in their traces and recorded profiles. The maximum absolute difference was 2.39e-12 and the maximum difference divided by max(1,absolute values) was 4.01e-13. This is a replay of the same implementation, not an independent discretization or external peer review.

A time-1 single-object pilot and a short time-reversal preflight preceded the freeze and are disclosed in [pilot.json](pilot.json). The supplementary manufactured-operator, stationary-profile, and initial-current checks were added after the completed physical evolutions. They and the figure generation ran in V8 on the retrieved data; Node commands are supplied for reproduction. They did not change the physical evolutions or their parameters.

## One unchanged classical model

Use tau=g v t, rho=g v r, z=g v z_physical, and X=v(u+i w)/sqrt(2), S=phi/v. The flat-space D=0 equations from the parent's action are

\[
\begin{aligned}
u_{\tau\tau}&=\nabla^2u-[S^2+b(u^2+w^2)]u,\\
w_{\tau\tau}&=\nabla^2w-[S^2+b(u^2+w^2)]w,\\
S_{\tau\tau}&=\nabla^2S-[a(S^2-1)+u^2+w^2]S,
\end{aligned}
\]

with cylindrical axisymmetric Laplacian

\[
\nabla^2=\partial_\rho^2+\rho^{-1}\partial_\rho+\partial_z^2.
\]

All physical evolutions use a=1 and b=0.1, equivalent to g_X=0.1, lambda_phi=0.01 and lambda_X=0.001. The reporter remains D=0. Gravity is not evolved. This is the same classical flat-space scalar sector used by the parent's weak-gravity benchmark, with no fitted pair force, external damping, or collision-specific coupling.

The energy and positive-charge convention are

\[
\begin{aligned}
\bar E={}&\int d^3x\left\{
\tfrac12(\dot u^2+\dot w^2+\dot S^2+
|\nabla u|^2+|\nabla w|^2+|\nabla S|^2)
+\tfrac a4(S^2-1)^2
+\tfrac12S^2(u^2+w^2)
+\tfrac b4(u^2+w^2)^2\right\},\\
\bar Q={}&\int d^3x\,(w\dot u-u\dot w).
\end{aligned}
\]

Physical quantities satisfy E=(v/g_X) Ebar, Q=Qbar/g_X^2, and physical length/time equal dimensionless length/time divided by g_X v. At g_X=0.1, multiply Qbar by 100 for the parent's physical charge convention. These are classical field charges, not exact quantum charge-eigenstate preparations.

Dimension, metric background, local field equations, initial charged resources and axisymmetry are supplied. No previous moving-word or M1 memory theorem is imported.

## Preparation and run matrix

The radial F,S profile is read from the parent's archived CSV (Git blob 0e221c32e604cf2ff89f0de8edc60ba4ef417ddd). Every tenth row is retained, giving 1,001 points on radius 0 through 50. Cubic Hermite interpolation uses the stored derivatives. Sampled differences against the native table are below 4.24e-9 for fields and 2.64e-7 for derivatives. Outside radius 50 the fields take vacuum values, matching the archived reporter convention.

A stationary object has u+i w=F(r) exp(-i Omega tau), with Omega=0.8000000019763136. For each moving object, the input is the Lorentz boost of that radial profile. The two initial centers are z=-12,+12 and velocities +0.2,-0.2. The complex fields are added; S=1+(S_left-1)+(S_right-1). Canonical velocities include both profile motion and phase rotation.

This is a supplied two-object initial state, not an exact stationary two-object solution or an autonomous assembly result. Interaction energy is included in the evolved initial Hamiltonian. Each grid samples the same continuum preparation; the discrete charge differs slightly through quadrature/discretization and approaches the continuum value on refinement.

At each of h=0.5, 1/3 and 0.25, the suite contains:

1. A static isolated object.
2. A boosted isolated object.
3. Equal-charge objects with relative phase 0.
4. Equal-charge objects with phase pi/2.
5. Equal-charge objects with phase pi-0.2.
6. Equal-charge objects with phase pi.
7. Equal-charge objects with phase pi+0.2.
8. Opposite-charge objects, with opposite rotation frequencies and the specified zero relative phase.

All evolve to tau=80. Timesteps are 0.1, 1/15 and 0.05 respectively. The standard domain has radial extent 48 and axial half-length 64. One additional medium-grid phase-zero run halves the timestep. Two additional medium-grid runs enlarge the domain to radius 56 and axial half-length 72, at phases zero and pi. Total: 27 distinct prescribed evolutions, plus the eight coarse runtime replays.

The phase-neighborhood runs were prescribed before collision outcomes were calculated. Exact phase pi preserves a central antisymmetry node; its behavior alone cannot establish a generic separation mechanism.

## What the completed runs show

![Finest-grid axial profiles and charge-sharing histories](RESULTS.svg)

The table reports the finest grid at tau=80. Positions are centroids of axial regions above five percent of the maximum axial F-squared line density. They describe field concentrations, not certified particle centers or asymptotic states.

| Prepared case | Endpoint axial structure | Approximate concentration positions | Charge in left / right half |
|---|---|---|---|
| Equal charge, phase 0 | One central connected density region | Centered near zero | 50% / 50% |
| Equal charge, phase pi/2 | Two separated outgoing density regions | -10.9, +5.9 | 33.44% / 66.56% |
| Equal charge, phase pi-0.2 | Two separated density regions | -9.5, +8.8 | 47.69% / 52.31% |
| Equal charge, phase pi | Two separated density regions | -9.1, +9.1 | 50% / 50% |
| Equal charge, phase pi+0.2 | Reflection of the pi-0.2 case | -8.8, +9.5 | 52.31% / 47.69% |
| Opposite charge | A connected, structured central distribution | No single-particle assignment | Net charge approximately zero; local positive/negative charge remains |

Changing the axial density threshold to two or ten percent gives the same component counts in all finest-grid cases. That is a useful diagnostic robustness check; it is not a proof of particle identity or nonlinear stability.

In the quarter-cycle case, the left-half Qbar changes from 399.7061 to 267.3092 while total Qbar remains 799.4166. The net transfer is about 132.397 Qbar, or 33.12% of the initial left-half charge. The outgoing threshold regions together contain about 97.0% of the total F-squared norm. At the five-percent threshold their charges are about 251.3 and 521.6 Qbar; remaining charge is in the surrounding fields. Half-space charges therefore should not be mistaken for exactly isolated object charges.

The charge shares are still changing near the endpoint: left-half Qbar is about 271.85 at tau=60, 276.35 at 70 and 267.31 at 80. The result is finite-time exchange and separation, not a settled asymptotic charge partition. The equal-phase central state also has substantial transverse extent and remains dynamical; a final stable fused soliton has not been established.

The opposite-charge case requires additional caution. It has the strongest discretization sensitivity of the suite. Net charge conservation does not identify an annihilation fraction, and a connected F-squared region does not prove one surviving object. No clean annihilation cross section, final state, or lifetime is claimed.

## Numerical controls and remaining error

The radial finite-volume cells are centered at rho_i=(i+1/2)h. Their weights are m_i=2 pi rho_i h^2. Radial flux coefficients are (i+1)/(i+1/2) and i/(i+1/2); the latter vanishes at the axis. Axial differences use the same spacing. Fixed vacuum ghost values terminate the outer grid.

The spatial force is the derivative of a matched weighted bond Hamiltonian. Interior bonds are counted once in the energy, with the same coefficients as the Laplacian. Boundary ghost bonds are included. Kick-drift-kick Verlet then gives a reversible symplectic discretization. Its global U(1) charge is preserved by both substeps to roundoff. Energy has a timestep error and is measured rather than assumed conserved.

| Check | Observed result |
|---|---|
| Vacuum control | Exactly unchanged in the recorded two-unit test |
| Forward/reverse control | Hosted replay: maximum field difference 3.34e-15; momentum difference 1.14e-14 |
| Global charge, all 27 cases | Maximum error divided by initial absolute charge 1.26e-13 |
| Matched energy, all 27 cases | Maximum drift 7.06e-4, or 0.0706%; largest in the coarse opposite-charge case |
| Fine quarter-cycle energy | Maximum drift 5.59e-6, or 0.000559% |
| Initial isolated Ebar, coarse to fine | Relative continuum-reference error decreases from 0.1262% to 0.0313% |
| Initial isolated Qbar, coarse to fine | Relative error decreases from 0.2711% to 0.06755% |
| Boosted isolated center at tau=80 | Finest result 3.9682 versus continuum kinematic value 4 |
| Manufactured three-dimensional operator | For F=r^2, initial u_tt=6-r^2 is reproduced to 2.01e-13 over the tested interior points |
| Stationary-profile acceleration defect | Relative X defect falls from 0.00574 to 0.00145 between h=0.5 and 0.25, consistent with second-order spatial error |
| Phase-neighborhood reflection | The two pi +/- 0.2 runs agree after axial reflection within the declared roundoff-level comparison |

The manufactured operator check uses a=b=0 only to isolate that control; none of the 27 physical evolutions changes the benchmark couplings.

Across the five recorded profile times, the largest normalized axial-profile L1 difference between medium and fine is about 0.99% for phase zero, 1.65% for pi/2, 1.12% for pi, and 2.09% for opposite charge. The medium-to-fine differences are smaller than the corresponding coarse-to-medium differences. These are convergence observations, not rigorous continuum error enclosures.

The opposite-charge transient has up to about 4.66% medium-to-fine difference in total F-squared norm and 0.0266 difference in central norm fraction. Its fine-scale time history needs further refinement before detailed channel claims.

Halving the medium-grid timestep in the phase-zero case changes the recorded normalized axial profiles by at most 0.184%. Enlarging the box changes the same measure by at most 1.95e-7 for phase zero and 6.95e-8 for pi. Boundary changes are much smaller than the tested spatial discretization changes over this window. They do not establish boundary independence at arbitrarily late times.

Axisymmetry permits three-dimensional head-on dynamics but excludes off-axis trajectories and non-axisymmetric instabilities. No initial-state distribution or statistical confidence interval was supplied. The sparse phase samples do not prove an open interval of stable behavior, a thermodynamic phase, or a populated cosmology.

## Why relative phase moves charge: an analytic bridge

The charge convention above has continuity equation

\[
\partial_\tau q+\nabla\cdot\mathbf j=0,\qquad
q=w\dot u-u\dot w,\qquad
\mathbf j=u\nabla w-w\nabla u.
\]

This is the ordinary Noether current of the supplied complex scalar action. It is not a new fundamental principle.

For the supplied symmetric initial two-object superposition, take half-separation d, opposite velocities +/-v, gamma=(1-v^2)^(-1/2), and right-minus-left phase theta. At the midpoint z=0 the two radial amplitudes agree. Their advective current contributions cancel, leaving

\[
j_z(\rho,0)=
-\frac{2\gamma^2 d}{r}\,F(r)F'(r)\sin\theta,
\qquad r=\sqrt{\rho^2+\gamma^2d^2}.
\]

Integrating over the midpoint plane, using F(infinity)=0, gives the initial relation

\[
\boxed{\left.\frac{d\bar Q_{\rm left}}{d\tau}\right|_{\tau=0}
=-2\pi\gamma^2d\,F(\gamma d)^2\sin\theta.}
\]

For this preparation, F(gamma d)=0.00203871371224 and the coefficient is 0.000326439248426. Its value follows from the same radial profile used for the isolated object, with no collision-specific fit.

The independently evaluated finite-volume midpoint current at phase pi/2 approaches this coefficient:

| h for the initial midpoint-current quadrature | Discrete coefficient | Relative difference from continuum expression |
|---|---:|---:|
| 0.5 | 0.000332700386 | 1.918% |
| 1/3 | 0.000329211862 | 0.8494% |
| 0.25 | 0.000327996858 | 0.4772% |
| 0.125 | 0.000326828173 | 0.1191% |
| 0.0625 | 0.000326536448 | 0.02978% |

The last two rows are initial-current quadratures, not additional time-80 evolutions. The phase-zero current vanishes to roundoff; the pi +/- 0.2 currents reverse sign as predicted. This relation explains the initial direction of transfer. It does not replace the full nonlinear evolution with a constant sine coupling throughout the collision. Profiles and phases change substantially during the encounter.

A complementary exact discrete statement is that each internal bond transfers equal and opposite charge to its two cells. For the bond across the midpoint plane,

\[
\frac{d\bar Q_{\rm left}}{d\tau}
=\sum_i 2\pi\rho_i\,(w_Lu_R-u_Lw_R).
\]

This is the flux used in the supplementary check. It is consistent with the weighted Hamiltonian and clarifies how local charge can change while the closed grid's total remains conserved.

## What this earns for the working theory

The earlier repository review identified a missing coexistence test: do objects retain distinguishable organization while interacting under the same microscopic law? This calculation supplies a concrete finite-time example for prepared classical scalar objects. Their relative phase affects whether the endpoint is central or separated and how conserved charge is redistributed. No new pair interaction or sector-specific rule was introduced.

It also shows why the lower energy of one large static lump does not settle collision outcomes. Dynamical state information, including phase and incoming motion, matters. Some tested channels remain separated over this window even though the compared fused endpoint can be energetically favorable. That is not a permanent barrier or a many-object stability theorem.

The broader closed-whole question remains open. Closure enforces common accounting of charge, energy and transferred information. It does not fix the incoming relative phase, select the scalar action, or derive the prepared objects. Ordinary interacting field theory already contains the mechanism calculated here.

The next scientifically distinct work is to establish post-encounter object stability and sensitivity to non-axisymmetric perturbations, then test formation from a stated initial ensemble if this model remains the chosen anchor. A reporter/record calculation would be separate: D=0 here, so no signal instrument or lasting acquired record was constructed. None of these results supplies quantum outcome selection, gravity beyond the adopted framework, or subjective consciousness.

## Reproduce and inspect

From this directory:

~~~bash
node run.cjs
node verify.cjs
node additional_controls.cjs
node plot.cjs
~~~

The main suite requires only Node's standard library. Running a single named case is also supported, for example:

~~~bash
node run.cjs fine_pair_phase_halfpi
~~~

The full verifier expects all 27 case files. Each records its configuration, grid, 81 one-time scalar diagnostic rows and five axial profile frames. These are the raw diagnostic outputs; full two-dimensional field snapshots at every step were not retained. The figure plots those recorded values directly. The source kernel and archived profile are sufficient to rerun the deterministic preparations.

The branch-scoped workflow leaves the ordinary builder branch and default branch unchanged. Historical source and verdict files are preserved, apart from the coordination entry describing this separate research package. No pull request, merge, external message, independent reviewer, or empirical result is claimed.

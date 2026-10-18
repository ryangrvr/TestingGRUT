# Adversarial challenge: Candidate 1 finite-time stability tests

18 October 2026. This audit holds the Candidate 1 action and physical couplings fixed. It introduces no replacement mechanism and inherits no result from the earlier moving-word architecture. Initial review concerns the proposed branch/fragment comparison and nonlinear three-dimensional evolution; the final output audit is recorded below when the results exist.

## Scope and rejection authority

These calculations may reject a selected preparation's claimed classical stability by showing converged breakup, persistent deformation growth, charged emission, collapse of its binding measure, or an energetically favorable allowed charge partition. A finite collection of surviving trajectories cannot prove infinite-time stability, nonlinear stability against every perturbation, a quantum lifetime, or the absence of nonradial modes not excited by the tests. A failed numerical conservation or boundary control invalidates that run's physical interpretation.

The strongest acceptable positive outcome is a stated finite-time response of the same-charge preparation to the explicitly listed perturbations, with quantified discretization and boundary errors, together with finite-charge energy comparisons over the sampled branch. Negative results must be retained.

## Dimensionless normalization

Set `A=sqrt(2) X/v`, `S=phi/v`, `tau=g_X v t`, and `xbar=g_X v x`. With `lambda_phi/g_X^2=1` and `lambda_X/g_X^2=0.1`,

\[
\ddot A=\bar\Delta A-[S^2+0.1|A|^2]A,
\qquad
\ddot S=\bar\Delta S-[(S^2-1)+|A|^2]S.
\]

The dimensionless potential, energy, and positive charge are

\[
\bar U=\tfrac14(S^2-1)^2+\tfrac12S^2|A|^2+0.025|A|^4,
\]

\[
\bar E=\int d^3\bar x\,[\tfrac12(|\dot A|^2+|\bar\nabla A|^2+
\dot S^2+|\bar\nabla S|^2)+\bar U],
\quad
\bar Q=-\int d^3\bar x\,\operatorname{Im}(A^*\dot A).
\]

Conversion back to physical units is `E=(v/g_X) Ebar` and `Q=Qbar/g_X^2`. The rotating profile is `A=F exp(-i Omega tau)`, giving positive `Qbar=Omega integral F^2`. In particular, `E/(m_X Q)=Ebar/Qbar` because `m_X=g_X v`. The evolution equations and the proposed normalization are consistent.

## The discretized Hamiltonian must match the evolution operator

For a real symmetric negative-semidefinite discrete Laplacian `L_h`, use the same operator in the force and the Hamiltonian. The gradient energy on a periodic grid is

\[
\bar E_{\rm grad,h}=-\frac{h^3}{2}\sum_n\left[
\operatorname{Re}(A_n^*(L_hA)_n)+S_n(L_hS)_n\right].
\]

Squaring a separately chosen centered first derivative does not in general give the Hamiltonian of the standard fourth-order centered Laplacian. That mismatch can create spurious reported energy drift or conceal an evolution error. FFT evolution should use its matching spectral `k^2` energy; the nonlinear collocation potential must use the same sample/quadrature convention. Aliasing and resolution still need continuum controls.

Kick-drift-kick leapfrog preserves the finite-grid U(1) angular momentum to floating-point accuracy when its force is derived from a U(1)-invariant discrete Hamiltonian. The reason is zero total U(1) torque during kicks and unchanged angular momentum during free drifts. Consistently paired staggered diagnostics can also represent this invariant. Define the actual diagnostic, rather than reconstructing velocities arbitrarily. Leapfrog is symplectic but does not conserve its original discrete energy exactly; bounded energy error and its time-step convergence must be demonstrated separately.

The massless three-dimensional fourth-order Laplacian has `max(-L_h)=16/h^2`; leapfrog therefore requires `dt<h/2` even before adding the potential Hessian. For a cubic FFT grid the corresponding upper bound is approximately `dt<2h/(pi sqrt(3))`. Include the largest positive coupled local-potential Hessian eigenvalue when choosing the actual step. These are stability limits, not accuracy targets. Refine dt at fixed h at least once so spatial and temporal error are not inferred from a single coupled refinement.

### Actual selected integrator: exact linear subflow with nonlinear kicks

The root subsequently selected Fourier collocation with Strang splitting, not the finite-difference/leapfrog option. Write `q=(Re A, Im A, S-1)` and use

\[
H_{\rm lin}=\tfrac12\int[p^2+|\nabla q|^2+q_1^2+q_2^2+2q_3^2],
\quad
V_{\rm nl}=\bar U-\tfrac12|A|^2-(S-1)^2.
\]

The exact Fourier Klein–Gordon subflow has frequencies `sqrt(k^2+1)` for A and `sqrt(k^2+2)` for `S-1`. The local nonlinear kicks are

\[
\Delta\dot A=-\Delta\tau[(S^2-1)+0.1|A|^2]A,
\]

\[
\Delta\dot S=-\Delta\tau[(S^2-1)S+|A|^2S-2(S-1)].
\]

These formulas are correct. Both subflows preserve the same collocation charge `h^3 sum(q_2 p_1-q_1 p_2)`. The earlier leapfrog linear-CFL bounds do not constrain this exact linear subflow. Nonlinear splitting accuracy, possible step resonances, and continuum/aliasing error still require the proposed dt and h refinements. The exact linear update does not make the full nonlinear finite-step trajectory exact or conserve its original H exactly.

The proposed matched-spacing domain pair is `L=64,N=96` and `L=85.333...,N=128`, both `h=2/3`, with separate dt halving and fixed-box spatial refinements. This is a useful comparison design. Since the target duration `T=40` exceeds half the smaller box, the finite-box qualification must remain explicit and be decided by actual domain comparisons.

## Same charge means the same initial conserved resource

The unperturbed continuum profile interpolated onto a Cartesian grid is not generally a discrete stationary state. Its initial Q and E must be measured with the same grid quadrature used throughout evolution, and the unperturbed run measures interpolation/relaxation radiation as a control.

The perturbed and control preparations must share the intended Q within an explicitly reported tolerance. An amplitude quadrupole with zero angular mean changes Q at second order, so first-order orthogonality is not enough. Adjust the angular momentum/momentum or amplitude explicitly and report the adjustment. The injected energy and changes in other conserved quantities are part of the perturbation, not a free resource.

A smooth phase perturbation `A -> A exp(i epsilon chi)`, `A_dot -> A_dot exp(i epsilon chi)` preserves Q pointwise while adding spatial-gradient energy. An ell=2 amplitude perturbation must have a regular radial factor proportional to `r^2` near the origin; `F(r) Y_20` with nonzero `F(0)` alone is not smooth Cartesian data. Prefer a smooth Cartesian polynomial times a localized radial envelope. Avoid unintentionally exciting a translation, total momentum, or angular momentum unless it is deliberately part of the test.

Changing the charge, physical couplings, vacuum, or solver branch to evade an instability changes the tested preparation. Numerical continuation, resolution changes, and explicit preparation changes can be allowed, but the original failure remains in the record and comparisons must identify the actual state tested.

### Pre-trajectory amendment: volume-preserving ellipsoidal strain

Before executing the three-dimensional trajectories, the root selected the constructor's calculable family `F_epsilon(x)=F(|M x|)`, `S_epsilon(x)=S(|M x|)`, with `M=diag(exp(epsilon),exp(-epsilon),1)` and epsilon 0.05 or 0.20. The proposed earlier amplitude/phase family was not executed and must not be presented as tested. A dated pre-trajectory amendment preserves this methodological history without changing the action or couplings.

Since `det M=1`, the continuum charge, rotational time energy, and potential energy are unchanged when the initial momentum is transformed consistently. Spherical symmetry makes the original integrated gradient tensor isotropic, so

\[
K_\epsilon=\tfrac13\operatorname{tr}(M M^T)K_0,
\qquad
\Delta\bar E=\tfrac43K_0\sinh^2\epsilon.
\]

This expression is correct and supplies an independent initialization check. Both scalar profiles must be deformed. A finite-grid charge renormalization of the initial angular momentum changes its time energy and must be reported separately from the exact continuum strain energy. By parity this preparation has zero initial net momentum. It is quadrupolar at first order in epsilon, with higher angular harmonics at finite epsilon; it is not a pure linear ell=2 eigenmode. Its positive energy cost establishes a restricted deformation-direction result, not positivity of the full fixed-Q Hessian.

## Boundary control

FFT imposes periodic images. Finite differences with fixed or periodic boundaries impose their own box dynamics. Initial exponential tails must be negligible at the faces, and the observation time must be compared with the distance from the perturbed/source region and measurement region to the faces, with a stated propagation margin. `T<L/2` without radii is not a sufficient argument.

The continuum principal cone has speed one in these coordinates, but FFT has a nonlocal numerical stencil and leapfrog changes numerical dispersion. A continuum light-cone statement alone does not certify absence of finite-box artifacts. A larger-box run at the same h and dt, comparing central observables over the stated interval, is the preferred direct control. If absorbing/damping layers are used, energy and charge losses must be accounted as boundary flux/reservoir losses, not mistaken for physical soliton decay.

## What to measure

Global E and Q conservation are necessary but insufficient: a completely dispersed state can conserve both. Measure charge and energy retained in stated core radii, emitted flux or outer-region content, moments/quadrupoles, peak and radial profiles, and shape changes against the unperturbed control. Report the measurement radius and sensitivity to its choice; a changing profile can cross an arbitrary radius without actual escape.

Translations and a global U(1) phase are neutral symmetries. For shape comparison one may align center and global phase, but report the raw displacement and phase as well. Alignment must not remove a growing quadrupole or hide transport of charge out of the core. Compare perturbation growth after removing the measured unperturbed-grid drift, with raw diagnostics still available.

Require a quantitative error scale from h, dt, and box refinements. A diagnostic change comparable to those errors is inconclusive. A small absolute conservation error is not enough if the claimed growth/energy difference is smaller still. A reliable instability should persist under refinement and exceed the unperturbed control and estimated numerical uncertainty.

## Branch and fragment-energy tests

For a smooth stationary branch, `dE/dQ=omega`. In dimensionless units `dEbar/dQbar=Omega`. A useful diagnostic is negative `dQ/domega`, but its sufficiency for multi-field perturbation stability needs spectral assumptions. No exact constrained Hessian spectrum has been supplied merely by calculating that slope.

`Ebar/Qbar<1` blocks complete dissolution into free X charge at leading order. Fission or emission requires the stronger comparisons

\[
\bar E(\bar Q)<\sum_i\bar E_{\min}(\bar Q_i)+\bar Q_{\rm free},
\qquad \sum_i\bar Q_i+\bar Q_{\rm free}=\bar Q.
\]

Fragment energies must refer to the same action, asymptotic vacuum, matching convention, and physical charge. Charge fractions should be solved directly or interpolated only with quantified errors; a favorable sign smaller than interpolation or profile-energy error is unresolved. A branch solver can miss a lower-energy disconnected or nonspherical state. Restrict the conclusion to the branch and partitions actually tested. The positive thin-wall surface term motivates resistance to large-fragment fission but does not replace this finite-Q test.

Quantum charge is integer. Continuous-charge classical comparisons approximate the large-charge limit; conversion `Q=Qbar/g_X^2` should be stated when discussing few-particle emission. Neutral radiation can lower excess energy without changing Q. A preparation whose energy exceeds the branch minimum may relax without demonstrating instability of that minimum.

## Gravity and quantum limitations

The proposed three-dimensional equations are nongravitating leading scalar equations. Their evolution cannot discover gravitational collapse or certify a gravitating branch. Keep the previously computed compactness as an external regime diagnostic, and do not extrapolate to arbitrary charge. Classical grid evolution also cannot determine a quantum lifetime, tunneling rate, quantum depletion, or Hartree backreaction. The known FLS quantum warning is relevant but its numerical instability window does not transfer to this positive-self-quartic action.

## Output audit

### Branch results inspected directly

`branch_scan.py`, `branch_protocol.json`, and the structured `branch_results.json` were read. The reference is

\[
\Omega_*=0.800000002,\quad \bar Q_*=399.4384624,\quad
\bar E_*=363.9229986,\quad \bar E_*/\bar Q_*=0.91108652.
\]

All 75 logged BVP attempts were accepted by the code's declared stationary-profile criteria. This does not imply every possible seed/branch has been exhausted. The energy integrals use the correct normalization and include the S and F gradients, rotational energy, and full potential. The fixed-charge fragment BVP includes an integral charge auxiliary variable and re-solves Omega; it does not substitute an unmatched fraction of the original profile.

The resolved sampled comparisons have positive energy costs relative to the reference:

| Comparison | Added dimensionless energy |
|---|---:|
| Complete free-carrier dissolution | 35.51546387 |
| Two equal found spherical fragments | 30.67384256 |
| Found 40:60 spherical fragments | 29.33901713 |
| Half-charge found lump plus free carriers | 33.09465321 |
| Three-quarter-charge found lump plus free carriers | 18.56408461 |

The fragment energy changes under the reported tighter BVP/larger-domain refinement are at most approximately `4.1e-7`, far below these resolved margins. Thus the sampled comparisons do not reveal an energetically favorable fragmentation channel.

**Unresolved comparisons remain unresolved.** Three, four, or five equal fragments, and partitions involving 20%, 25%, or one-third charge lumps, were not solved: those target charges lie below the minimum of the found branch. Their absence from that branch is not proof that no other state, fragment type, or disconnected branch exists. The code appropriately marks these comparisons unresolved.

The turning-point minimizer explicitly reports `optimizer_success=false` because its maximum function-call budget was reached. Its confirmed stationary BVP near `Omega=0.95980423`, `Qbar=148.73404` is a useful approximate location of the found minimum. It is not a converged optimizer certificate or rigorous exact existence boundary. The separately located binding crossing is converged near `Omega=0.90654404`, `Qbar=176.60439`. Neither threshold is a new quantum stability boundary.

A subsequent, explicitly retained refinement was also inspected in `branch_thresholds_refined.json`. It solves the sensitivity condition `dQbar/dOmega=0` and now reports a converged root at `Omega=0.95980474869`. The resulting minimum is `Qbar=148.734042107`, with `Qbar` and `Ebar` differences of approximately `1.7e-9` between domains 80 and 100. The binding crossing is refined to `Omega=0.90654404078`, `Qbar=176.604388320`, with domain differences below `9e-10`. The first-pass optimizer failure remains historical evidence; the location is no longer an unresolved numerical optimization objection after this refinement. It still describes the found branch rather than every possible solution family.

The branch identity `E-Omega Q=(2/3)K` and `dE/dQ=Omega` imply decreasing E/Q along each smooth branch as Q increases. The code's extension from a found continuous energy envelope to all positive-charge partitions drawn from that envelope is logically conditional. The premise that the envelope includes every admissible lower-charge component is not established by a finite nodeless spherical scan. No global fixed-charge minimizer or all-partition theorem follows.

The completed second-pass data were subsequently reread: there are now 112 accepted logged profile BVP attempts, with the first-pass records retained. Additional exact-charge targets and both found-branch combinations remain energetically unfavorable. The smallest resolved sampled cost across the expanded tables is `+15.10481616`, for an 80%-charge lower-branch lump plus 20% free charge. The `+18.56408461` value in the first-pass table above is not the minimum of this enlarged set. An arbitrarily small emitted classical charge has an arbitrarily small energy cost; the finite sampled minimum is not a uniform fission/emission gap for every partition.

The radial fixed-Q Hessian calculation in the added branch code is consistent with the transformation `u=sqrt(4 pi) r delta F`, the Dirichlet endpoint conditions, the coupled local potentials, and the quadrature-weighted rank-one charge correction. Its lowest eigenvalue changes from `+0.0304422` to `+0.0305207` as N increases 800 to 1200 at Omega 0.95, and from `-0.0274040` to `-0.0273358` at Omega 0.97. This supplies a resolved radial energetic negative control on opposite sides of the turning point. It is a Hessian result for the stated amplitude class, not a complete rotating Bogoliubov spectrum or a nonlinear three-dimensional trajectory of the high-frequency branch. A negative eigenvalue here must not be converted directly into an asserted exponential growth rate. No full nonradial constrained spectrum has been computed for the reference.

### Independent finite-volume method check

`independent_notes.md` and `independent_refinement_results.json` were inspected. The eliminated-charge term `Qbar^2/(2I)` is correct. With volume-weighted field variables the displayed local Hessian plus `4 Omega^2 |u><u|/I` rank-one term is the correct full fixed-charge Hessian; omitting it would be wrong. The outer half-cell boundary coefficient and origin zero-face-area treatment match the described cell-centered radial functional.

All twelve original bounded L-BFGS-B runs stopped on relative-energy ftol rather than the requested projected-gradient tolerance. That limitation is retained. The later frozen Newton check of the three fine Gaussian-origin states reaches full gradient infinity norms from `4.37e-12` to `9.14e-12`; the tiny projection corrections do not account for this improvement. This supports actual discrete stationarity, rather than using close energies alone as a substitute for a stationary solution.

At charge fractions 1, 0.75, and 0.5, halving the cell spacing reduces the discrepancy from the BVP energies by approximately four. The fine-grid equal-halves gap is `+30.67374232`, consistent with the BVP `+30.67384256` and far from zero. The assumed leading-h-squared extrapolations agree with BVP values at approximately `1e-7`. This is useful agreement of methods/discretizations for these radial stationary energies. Both finite-volume grids use the same radius 40; they add no independent domain refinement. Two grids do not supply a rigorous continuum error enclosure.

No independent contradiction was found. The word independent here denotes a separate implementation and method within the same team. It is neither external validation nor empirical evidence. Newton stationarity does not prove Hessian positivity or global minimization, and the bounded radial nonnegative-profile class does not exhaust all three-dimensional field states.

### Nonlinear three-dimensional outputs

All seven declared outputs were subsequently inspected directly, including the full diagnostic histories and saved final slices. The tested interval is `0<=tau<=40`, sampled at integer tau. The fixed charge is the original `Qbar=399.4384624223012`; the on-grid frequency adjustments are explicit. No physical coupling was changed between runs.

`evolve.py` was read before interpreting its outputs. Its actual split, exact Fourier oscillator update, Q sign, on-grid angular-frequency normalization, and matching Parseval global Hamiltonian are consistent. Local regional gradient energy uses separately reconstructed spectral derivatives. On an even collocation grid Nyquist derivative conventions can make this differ slightly from the global k-squared Hamiltonian; the recorded `spectral_local_gradient_total_difference` must be checked, particularly on the coarse grid, before interpreting regional energies as resolved continuum densities.

The predeclared relative energy-drift gate `2e-3` is a broad numerical failure gate, not a sufficient resolution argument by itself: it would allow roughly 0.73 energy units, greater than the mild strain's approximately 0.222 energy injection. The actual energy drift, time-step differences, and unperturbed motion must be compared to the perturbation energy and response scales. The strong strain injects approximately 3.60 energy units. Satisfying a broad validity flag does not replace these quantitative comparisons.

The actual completed results are:

| Run | N, box side, dt | Minimum Q fraction inside r=8 | Maximum relative energy drift | Maximum phase-aligned A profile error |
|---|---|---:|---:|---:|
| Unstrained control | 96, 64, 0.1 | 0.997716 | 1.3201e-5 | 0.002207 |
| Strain epsilon=0.05 | 96, 64, 0.1 | 0.996962 | 1.6415e-5 | 0.049683 |
| Strain epsilon=0.20 | 96, 64, 0.1 | 0.984644 | 8.8389e-5 | 0.197242 |
| Time-step control | 96, 64, 0.05 | 0.984807 | 2.0168e-5 | 0.197242 |
| Coarse grid | 64, 64, 0.1 | 0.984570 | 8.8389e-5 | 0.197243 |
| Fine grid | 128, 64, 0.1 | 0.984749 | 8.8390e-5 | 0.197242 |
| Larger box, matched h | 128, 85.333..., 0.1 | 0.984644 | 8.8389e-5 | 0.197242 |

The maximum relative charge drift across every run is `6.12e-13`. Global conservation does not substitute for the regional results, but it excludes an unaccounted large source or sink in these trajectories. The largest relative discrepancy between the reconstructed local-gradient total and the global k-squared Hamiltonian is `6.63e-8` on the coarse grid; it is negligible at the quoted regional-energy precision.

The continuum analytic preparation energy is independently reproduced: the unstrained N96 difference is about `-1.49e-8`, the strong N96 difference is `+5.23e-7`, and the strong N128 difference is `+7.31e-10`. The coarse-grid difference is `-5.83e-4`, still small compared with the strong injected energy. This checks the intended strain and finite-grid normalization rather than inferring them from the desired final outcome.

Halving dt reduces the maximum energy error by a factor approximately 4.38, consistent with second-order splitting. The largest absolute energy error in the main strong case is approximately `0.0325`, about 0.9% of its injected strain energy. The mild run's approximately `0.00598` error is about 2.7% of its injected energy. The unstrained control has a real numerical motion floor: maximum A profile error `0.00221`, with final organization-profile error approximately `0.00591`. Deviations below that scale must not be interpreted as a demonstrated physical effect.

For the strong strain, the max sampled dt-halving differences are `1.62e-4` in the r=8 charge fraction, `3.34e-4` in A profile-error norm, `0.00115` in quadrupole norm, and `0.00370` in organization-profile error. Central S differs by as much as `0.00242`; numerical precision claims about the central depleted field must respect that scale. The N96-to-N128 spatial comparison changes the A profile-error history by at most `1.28e-8`. Hard-sphere mask quadrature changes the r=8 fraction by approximately `1.35e-4` between those grids; radius and histogram diagnostics should not be quoted as if they were exact continuum integrals.

The larger-box comparison is particularly informative. At matched spacing, the maximum difference is `1.21e-9` in core Q fraction at r=8 and `6.54e-8` in central S. The global A profile-error history differs by up to `1.21e-5`. An additional direct comparison of the saved final central slices at exactly matching physical coordinates gives differences at r<=12 no larger than `7.65e-8` for either A component or S. Over the full common square, however, outer-wave differences reach approximately `0.0024`. Thus perfect causal isolation is false as a numerical claim: the box affects outer content. The measured core response is nevertheless insensitive to this box enlargement over the tested interval.

**Direct field comparison is distinct from agreement of scalar diagnostics.** The small difference between two profile-error norms is not a norm of their solution difference. To check the latter, `compare_slices.py` spectrally resamples the final N128, L64, z=0 fields onto the N96 grid along both periodic axes and writes `convergence_slices.json`. At r<=12 the raw maximum differences are `3.58322e-6` in complex A and `5.71514e-6` in S. The relative plane L2 differences are `1.18920e-6` for A and `3.90151e-6` for S normalized to the organization profile S-1. The fitted global A phase is only `1.56e-9` radians; removing it changes these values negligibly.

Resampling truncation is recorded rather than hidden. Fine-to-reference-to-fine reconstruction has whole-plane maximum differences `2.71e-6` for A and `4.48e-6` for S, with relative plane L2 remainders approximately `1.06e-6` and `3.37e-6`. Thus the sampled spatial fields agree at a few parts in a million in this plane; the scalar diagnostic differences alone would have understated the direct field difference.

The same script repeats the exact-coordinate domain comparison. Inside r<=12 the raw maximum complex-A difference is `8.07e-8` and the S difference is `7.57e-8`; over the common full plane the complex-A difference reaches `0.00244`. Raw results are retained along with optional global-phase alignment, whose whole-plane fit can slightly worsen the tiny core difference because the outer fields differ. These are final two-dimensional slice comparisons, not a full three-dimensional solution norm or a rigorous continuum error bound. The time step is unchanged in the spatial comparison. Running the executable script regenerated the preexisting JSON byte for byte; no additional physical trajectory was run.

No sampled run crosses the predeclared 90% core-retention failure gate. The A profile-error and quadrupole norms in both strained runs remain no larger than their initial strain values over the sampled interval. The strong main run ends with A profile error `0.1681` and quadrupole norm `0.2152`, compared with initial values `0.1972` and `0.2648`; it has not returned to an exactly spherical stationary state. Its r=8 charge fraction declines from `0.996267` initially to `0.984644` at the final sample, and its r=12 fraction ends near `0.993032`. This is outward redistribution, not unchanged core occupation. It does not by itself prove permanent escape to infinity, and the decline is still present at the end of the window.

### Final adversarial verdict

**The selected reference survives this bounded challenge.** The resolved sampled fragmentation costs are positive; an independent fixed-charge discretization reproduces the radial energies; the radial Hessian check detects an energetic negative direction on the other side of the turning point; and all seven declared finite-time trajectories retain the scalar core with controlled conservation, time-step, spatial, and box comparisons. No observed result here requires rejecting the Candidate 1 action or this particular preparation over tau<=40.

The following stronger assertions remain rejected or unestablished:

- Global minimum at the reference charge, all possible fragments, all perturbation directions, or arbitrary initial data.
- Infinite-time nonlinear stability, complete relaxation, no outward redistribution, or a quantitative classical/quantum lifetime.
- A nonlinear Einstein evolution or exclusion of gravitational collapse at other charges.
- Quantum depletion, tunneling, Hartree stability, or an exact interacting coherent-state solution.
- Robustness of the original stationary reporter phase under the oscillating ellipsoidal target, repeated nondestructive readout, or memory fidelity. Those observables were not evolved here.
- A constructed laboratory source/detector or empirical relevance of the illustrative scalar sector.

The justified positive conclusion is finite-time persistence under the two specified, same-charge ellipsoidal strains of the unchanged leading scalar law. That conclusion is stronger than stationary-profile existence and substantially narrower than an unrestricted stability theorem.

# Independent expanding-charge cross-check, version 1

2026-10-21. Base: Codex `41e7dadc1afeb0eaf88e646aa91cf1f456e63f9d`.
Reviewed: Claude `c7b8480a49e3f0b8737e800396a30f33f6337981`, report, issue #13 reply,
charged-exit and bias notes, and selected separate-code verifier records. Their outcomes
were known before this protocol. This is a numerical cross-check, not a blinded
observation, an outside replication, or an empirical validation.

## Question and scope

Does the unchanged classical Candidate 1 action, with D=0, produce both vacuum signs
after dilution of a charge-held restored state in independently implemented 3D
evolution? Which finite-time conclusions survive seeds, box size, grid, timestep,
and the definition of a genuine vacuum region? A selected positive-vacuum, subcritical
charge preparation is a negative control. No new action, source, coupling, or
successful cosmological history is adopted.

This calculation can check a necessary link in a conditional wall constraint.
It cannot fix v, the primordial charge yield, the producer or the history. It cannot
establish persistence over cosmological times, calibrate a charged network's asymptotic
scaling law, or validate the extrapolated bias threshold near 10^-23. Equal marginal
probabilities alone will not be taken as proof of spatial percolation.

## Equations and units

Use ell = 1/(sqrt(lambda_phi) v), v=1 TeV, lambda_phi=g_X^2=0.01,
lambda_X/lambda_phi=0.1. Write X/v=(u+iw)/sqrt(2), S=phi/v.
Each real field has canonical kinetic normalization. Radiation expansion is externally
prescribed: a(tau)=tau/tau_i, tau_i=40, a(initial)=1. No gravitational backreaction,
thermal bath, D fluctuations, or quantum noise is included.

For chi=a(u,w,S), P=chi', the potential is

    V_tau = (chi_S^2-a^2)^2/4 + chi_S^2(chi_u^2+chi_w^2)/2
            + 0.025(chi_u^2+chi_w^2)^2.

Evolve chi''=laplacian(chi)-dV_tau/dchi with the original Cartesian solver's
seven-point periodic stencil and an independently written time-dependent KDK update.
No Claude solver code is imported or copied. The mathematical conformal-variable
method is shared and disclosed. Canonical charge uses the previous Codex convention:

    Q_com = integral(chi_w P_u-chi_u P_w) d^3x.
    n_physical = Q_com/(a^3 V_com).

The exact semidiscrete conformal-energy source is

    dH_conf/dtau = -a a' integral(chi_S^2-a^2) d^3x.

The spatially averaged physical density and pressure, including ALL simulated fields,
are a^-4 times K+G+V_tau and K-G/3-V_tau, respectively, where
K=mean(sum(P-a'/a chi)^2)/2 and G is the forward-link gradient energy. These are not
a unique split between walls and a receiver. No energy is discarded to manufacture
a wall component.

## Preparations

Restored: nbar_i=kappa sqrt(0.1), kappa=2 (main), kappa=4 (one additional control).
At S=0, F^2=u^2+w^2=kappa^(2/3), omega=sqrt(0.1) F. The homogeneous conformal
rotating orbit is exact. Add a declared finite, band-limited classical perturbation:
S rms 0.01; relative u/w perturbations 0.001; physical S velocity zero initially.
Fourier coefficients and phases are deterministic Gaussian draws; their complete
table is saved. Normalize the total initial charge to the stated value.

Cold control: kappa=0.9, positive broken branch, solve
q=F^2 sqrt(1-0.9 F^2) for its smaller positive root, S=sqrt(1-F^2).
Apply the same perturbation family. This preparation is an input, not charge production.

The initial perturbations are neither vacuum fluctuations nor a thermal distribution.
Grid refinements at fixed box and seed evaluate the same Fourier polynomial. Volume
comparisons hold correlation scale and rms fixed but are statistical, not nested copies.

## Frozen run list

All boxes are periodic cubes; L is the full comoving side. End at a=2.2 (tau=88).
Use seeds 101, 202, 303. The production suite is deliberately bounded:

| Group | N | dx | dt | seeds | kappa |
|---|---:|---:|---:|---|---:|
| main | 192 | 0.4 | 0.08 | 101,202,303 | 2 |
| larger box | 288 | 0.4 | 0.08 | 101,202,303 | 2 |
| finer grid, same L=76.8 | 256 | 0.3 | 0.06 | 101 | 2 |
| finest grid, same L=76.8 | 384 | 0.2 | 0.04 | 101 | 2 |
| half timestep | 192 | 0.4 | 0.04 | 101 | 2 |
| higher charge | 192 | 0.4 | 0.08 | 101 | 4 |
| cold control | 192 | 0.4 | 0.08 | 101 | 0.9 |

Before these runs, validate charge normalization, the homogeneous rotating orbit,
the linear crossing equation against a separate ODE, direct physical-field versus
conformal evolution, dt^2 convergence, exact Z2 mirror symmetry, periodic wrapping
on synthetic domains, and planar wall area/energy/velocity estimators. Small validation
runs may be refined when a defect is found; record the changes before production.

If cost exceeds 20 minutes of production wall time, stop launching jobs after the
current job and report the unexecuted rows. Do not replace failed runs or choose
seeds after inspecting outcomes. A resource interruption is not physical failure.

## Diagnostics and decision rules

Retain full-precision scalar time series, exact run configurations, mode tables,
source/environment hashes, selected field slices and reproduction commands. Formation:
each of S>0.5 and S<-0.5 occupies at least 10% of the volume at an output. Quote
the finite output-time bracket. Record fractions also at vacuum thresholds 0.4 and 0.6.

Compute six-neighbor connected components and noncontractible winding on the periodic
torus for both genuine-vacuum masks. Spanning a coordinate cut is not winding.
Record the largest component, winding axes, raw sign interfaces, and a basin interface
estimator that extends labels from genuine vacua through cores. Internal zero crossings
inside restored charged cores will not automatically be counted as walls. Topology
and interface output spacing is 4 conformal time units, plus the endpoint.

Report the dimensionless charge drift and the Hamiltonian source residual, with a
dt-halving control. Charge drift exceeding 1e-9 or energy-source residual exceeding
1% of max(initial H, current H, accumulated absolute work) fails the numerical gate.

Claim finite-window persistence only for the stated sampled times and seeds. Claim
resolution convergence only if successive refinements bring formation a within 0.05
and late scalar/area diagnostics within 10%, with a smaller change at finer resolution.
If topology changes with box or seed, report that dependence. No cosmological lifetime
or continuum percolation probability follows from a passed numerical gate.

## Interpretive controls

The one-speed Lorentz factor at rms speed 0.456 is NOT an upper bound on a distribution
of wall velocities. Distinguish area weights from energy weights. The proxy used in
the previous packet is not a calibrated energy-weighted velocity. Likewise an
instantaneous scalar pressure need not satisfy w<=1/3; a cosmological use requires
coarse-graining and a justified state family. A Delta N_eff-equivalent energy ratio
does not by itself constitute a Planck likelihood for arbitrary wall/matter stress.

Any defect, changed estimator or changed numerical protocol is a dated amendment.
The empirical-input problem remains open regardless of these run outcomes.

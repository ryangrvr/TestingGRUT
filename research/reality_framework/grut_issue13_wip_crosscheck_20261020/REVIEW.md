# Preliminary cross-check of Claude WIP 559975f

This reviews the published **unverified** modules at
`559975fad327672f70221e62f6c11c913f7d30b6`, in
`research/reality_framework/grut_issue13_wall_history_20261020/modules/`.
It does not incorporate the pending eight-checker workflow or charged-exit
simulation. No verdict on those unavailable results is implied. The seven
modules have not been independently re-executed here.

The substantive connection is sound at leading classical order: the same
interaction that depletes phi inside bound charged cores can stabilize
phi=0 throughout a sufficiently dense homogeneous charged preparation.
What happens as that preparation dilutes is a cosmological question, and
the previously tested cold fragments do not settle it. This is a useful
microphysics-to-history connection, not empirical validation of new fields.

One additional controlled piece is now explicit: the neutral planar wall
supports a localized X mode below the bulk threshold. This gives an
analytic matching check for the charged-exit work, described below. It does
not establish charge capture, an occupied-wall equilibrium or its lifetime.

## Findings requiring attention before synthesis

| Item | What the WIP actually establishes | Required qualification or calculation |
|---|---|---|
| Higgs portals remain exactly zero | Matter diagrams preserve zero portals when the sectors share no nongravitational vertex. hidden_heating explicitly estimates graviton-induced effects separately. | Retain **nongravitational RG** in every headline. An exact statement about the full gravity EFT does not follow. No new graviton-loop result is claimed by this review. |
| Inflation and gravitational-heating thresholds | Conditional stochastic/mode/rate calculations, with initial state, duration, reheating model and criteria supplied. | The 0.436v value concerns an equilibrium marginal, not all finite inflationary histories. The WIP itself gives different finite-duration and rare-domain criteria. Four quoted preparation conditions are not a proved necessary-and-sufficient classification of every possible history. |
| Second-order re-breaking | The approximately neutral, chemically equilibrated hidden-bath calculation gives a continuous transition at its approximation order. | Do not carry that conclusion into cold finite-charge dilution. The charge module itself finds a first-order structure, metastability and inhomogeneous coexistence. |
| Supply below 198 r^(-1/3) TeV avoids the charge restoration route | 198.20 TeV is the cold phi=0 **local spinodal** at r=1, v=1 TeV. | The global homogeneous switch, disappearance of the broken homogeneous minimum and coexistence endpoint are different thresholds (below). Initial state and charge localization matter. |
| Equal plus/minus probabilities imply a network | Z2-invariant ensembles have equal one-region sign marginals. The charge module correctly notes that this does not alone establish domains. | Match spatial sign correlations and the actual preparation to the lattice ensemble. Percolation thresholds for independent sites do not automatically apply to correlated inflationary or charged states. |
| One infinite wall carries ~99% of area; A=.78–.88 | Periodic finite-box, X=D=0 simulations have a box-spanning cluster. The most useful 3D PRS value is ~.777; the .884 physical-equation estimate is already box affected and has no box-safe scaling window after formation. | Use **box-spanning**. Keep PRS, physical equations, finite observation time, area estimator and energy estimator distinct. Do not label these charged-exit results. |
| Backreaction and T_dom=6.367 keV | cosmo_observable imposes rho_w=A sigma/t and integrates the Friedmann constraint, with a phenomenological post-domination continuation. | A wall-stress/energy-transfer matching remains to be supplied. Effective dilution pressure differs from microscopic wall pressure when energy is lost. See the deciding calculation below. No corrected physical T_dom is claimed here. |
| Restoring histories are all excluded | A specified persistent macroscopic network at the benchmark tension is strongly problematic under the declared cosmological continuation. | Establish the actual formation, loss products, persistence, vacuum/geometry matching and observable assumptions for the stated history. No universal exclusion or field identification follows from the unfinished synthesis. |

## Deciding cross-module check: wall pressure and energy losses

The WIP correctly derives the isotropic thin-wall relation

\[
p_w/\rho_w=\langle v^2\rangle_E-2/3.
\]

Here the average is weighted by wall energy. For a single speed, a planar
wall has energy per area gamma sigma. Intrinsic pressure, physical area
and the gamma factor must be matched to the lattice estimators. In
`network_sim/code/wallsim.py`, A is computed from weighted link **area**,
whereas the cosmology uses an **energy-density** A. A moving-planar-wall
control and an energy-weighted velocity calibration would establish that
matching. The one-speed gamma at the quoted v=.456113 is 1.12369; that is
a control value, not a demonstrated 12.4% correction to the whole network.

If we provisionally identify that velocity with the appropriate energy
average, it gives w_w=-.458627, whereas a separately conserved component
with rho_w proportional to 1/t in a radiation era would have w_eff=-1/3.
These can both describe a network **only if losses are accounted for**:

\[
\dot\rho_w+H(1+3\langle v^2\rangle_E)\rho_w=-Q,
\qquad
Q=\left[1-Ht(1+3\langle v^2\rangle_E)\right]\rho_w/t.
\]

For Ht=1/2 and the quoted velocity,

\[
Q=.1879412\,\rho_w/t=.3758824\,H\rho_w.
\]

The wall-only scaling component is therefore not isolated. In the complete
closed theory this energy goes somewhere: scalar particles, closed wall
fragments, gravitational waves or another explicitly modeled component.
The radiation/matter receiving it must have a corresponding positive source
in its evolution equation. Altering a pressure to reproduce dilution is a
phenomenological prescription, not by itself a derivation of the complete
stress tensor and emitted energy. The WIP cosmological core currently
includes ordinary radiation/matter and the prescribed wall density, but
no independently evolved receiving component for this scaling loss.

The one-scale literature explicitly includes loss terms; the background
scalar action also requires total covariant stress conservation. This check
does not say that the numerical Friedmann integration is arithmetically
wrong, that emission will rescue the benchmark, or that its large mismatch
must disappear. It identifies the matching needed before calling that
integration the complete backreacting history of the stated wall network.

The supplied independent constant-g ODE illustrates the effect. It keeps
rho_w=A sigma/t and the quoted constant velocity, evolves emitted matter
or radiation with the required Q, and stops at wall/visible-radiation
equality. The emitted-to-wall ratios there are .20566 (massive products)
and .08690 (radiation products). These are **toy consistency diagnostics**;
they neither update the WIP's numerical 6.367 keV value nor predict a
Candidate-1 population. True emitted spectra, charge, velocities, g_* and
network evolution require calculation. All energy transfers should be
matched before coupling the network and relic modules.

Similarly, no SM portal means no portal-mediated SM drag; it does not
automatically eliminate hidden-particle drag. Stable relics and emitted
quanta need not be Boltzmann suppressed at the **visible** temperature.
Their distributions and wall-scattering response must be calculated or
bounded. The special linear reflectionless result below is a useful
controlled case, but is not a proof for an interacting curved network.

Primary comparison: [Avelino, Martins & Oliveira, hep-ph/0507272v1](https://arxiv.org/html/hep-ph/0507272v1),
section III, especially equations 15–18. This review independently derives
the pressure/source identity and uses that paper for context. It does not
adopt its phenomenological fitted parameters.

## Independent cold-charge check

In the same conventions as the WIP,

\[
e(q,S,x)=q^2/(2x)+(S^2-1)^2/4+aS^2x/2+bx^2/4,
\quad a=1,\ b=.1,
\]

where q=n/(sqrt(lambda_phi)v^3), x=2|X|^2/v^2 and S=phi/v.
Directly eliminating S gives the two branches in the charge module.
An independent root solve reproduces:

| Threshold | n/v^3 | T_visible at r=1, v=1 TeV |
|---|---:|---:|
| Restored-phase local spinodal | .0316227766 | 198.198 TeV |
| Global homogeneous switch | .0401026192 | 214.531 TeV |
| Broken homogeneous radial spinodal | .0427666866 | 219.180 TeV |
| Dense endpoint of inhomogeneous coexistence | .1778279410 | 352.451 TeV |

The cosmological conversion retains the WIP's quoted baryon-number
calibration (Y_B=8.67400003e-11) and visible g_*s=106.75. It is not an
independent calibration of either r or v. Thresholds scale as v r^(-1/3)
at fixed dimensionless couplings and stated visible-entropy assumptions.

Near S=0 at fixed q, locally eliminating x gives the S^4 coefficient

\[
1/4-a^2/(12b)=-7/12.
\]

Independent finite differences reproduce it. It is negative, so the cold
finite-charge transition cannot be summarized as the neutral thermal
second-order transition. The charge module's distinct metastability and
coexistence analysis should survive into the synthesis. The local spinodal
is also not a proof that nucleation, fragmentation or sign-memory erasure
is absent at earlier times. In particular, applying a vacuum O(4) thin-wall
bounce at finite conserved charge needs justification; its estimate does
not close that dynamics question near the spinodal.

## Conditions for interpreting the pending charged-exit run

First, there is a useful exactly solvable planar control. With
phi=v tanh(z/delta), delta=sqrt(2/lambda_phi)/v, X=D=0, define u=z/delta.
The transverse fluctuation operators, after subtracting their asymptotic
mass, are

\[
L_\ell=-\partial_u^2-\ell(\ell+1)\operatorname{sech}^2u,
\quad \ell_X=1,\quad \ell_\phi=2,\quad
\ell_D=\tfrac12\left[\sqrt{1+8h_D^2/\lambda_\phi}-1\right]=.01961524.
\]

The X value follows from 2g_X squared/lambda_phi=2 at the exact tree
benchmark. The factorization
L_ell=(−partial_u+ell tanh u)(partial_u+ell tanh u)−ell squared
gives a localized X wavefunction proportional to sech u, with unshifted
eigenvalue -1. Restoring the positive bulk mass term gives

\[
\omega^2=k_\parallel^2+g_X^2v^2-\delta^{-2},
\quad m_{X,\rm surface}=m_X/\sqrt2=70.710678\ {\rm GeV}.
\]

This is a normal mode on an ideal uncharged planar wall. The negative
eigenvalue of the **subtracted** operator is a bound level, not a tachyon:
the physical squared rest frequency is positive. Independent finite-box
eigenvalues at du=.04 and .02 converge to -1. For phi, the two bound levels
are -4 and -1 before adding its bulk mass term 4: the translation mode
has zero frequency, and the shape mode has mass sqrt(3)/delta=122.47449
GeV. These are wall modes, not the old finite-composite transverse norm.

For integer ell, the same ladder operators map scattering states to free
plane waves without a reflected component. Direct complex-wave ODE
integration checks reflectionless linear phi and X scattering; flux errors
are below 8e-11 on the tested momenta. D is a noninteger control: its
reflection probabilities at p_normal delta=.05,.2,1 are .132284,.008366
and 2.8435e-5, matching the analytic potential formula. These are **linear,
static, planar, leading classical** results at the exact couplings.

Implications: the wall can support a low-amplitude charged X excitation,
and its ideal linear X scattering is transparent. A static linear
scattering solution does not populate the bound state: capture must dispose
of energy and momentum through additional dynamics. Running, wall curvature,
finite charge, nonlinear scattering and recoil can change this ideal result.
No zero-friction theorem, capture rate, stable finite charged wall,
cosmological lifetime or ordinary detector signal has been derived.

For the pending expanding charged simulation:

- Report the physical preparation, initial mean charge, spatial charge and
  noise distributions, plus/minus bias, field momenta and any residual
  coherent phi motion. A vanishing mean is not a specification of sign
  correlations or of quantum fluctuations.
- Distinguish initially homogeneous metastable charge from phase-separated
  charge, and quantify unresolved charge/energy outside identified cores.
- Record a^3 Q conservation and the full hidden stress-energy balance with
  FRW expansion. Conservation of charge alone does not verify the energy
  transfer or the geometry of a wall network.
- Count resolved interfaces between near-vacuum regions, with charged-core
  sign oscillations excluded or separately classified. Reuse the plateau
  and charged-field controls from the Codex packet, not sign(phi) alone.
- Match the quench/expansion rate. With the original r=1, v=1 TeV conversion,
  H/m_X at the 198 TeV local spinodal is **5.519e-10**. A numerically sped-up
  expansion can establish a mechanism in another preparation, but needs a
  controlled extrapolation or justified rescaling before predicting the
  original history. The microphysical wall width, growing correlation
  length, horizon, lattice size and time step all need to be identified.
- Distinguish finite bubbles and box-spanning walls; test resolution,
  duration, volume, noise spectrum/realization and coherent sign memory.
  Persistence over a lattice window is not persistence to 6 keV.
- Match the resulting network's area, energy, pressure, loss spectra and
  charge to the cosmology/relic evolution. The existing X=D=0 PRS network
  does not supply this matching automatically.

An explicit counterexample in check.py has equal plus/minus probability
at each of four cells in two Z2-invariant ensembles. The fully correlated
ensemble has zero sign interfaces; independent signs have neighbor
disagreement probability 1/2. This illustrates the missing joint-state
information only; it does not assert either ensemble is the actual thermal
or primordial state.

## Status and reproduction

Run from repository root:

```bash
python research/reality_framework/grut_issue13_wip_crosscheck_20261020/check.py
```

Only numpy and scipy are required. The script imports no Claude calculation
code and runs in a few seconds. RESULTS.json contains all outputs, source
hash, environment and 33 passing algebra/numerical checks. PROVENANCE.json
identifies the WIP source files and the exact read/computation scope.

No new lattice campaign, model extension, parameter fit, observational
comparison or claimed empirical validation was performed. All reviewed
headlines were available before these calculations, which are retrospective.
All original frozen packets and Claude's working branch are unchanged.
These checks do not substitute for the seven pending module verifications.
No new scientific numerical check failed; no failed evolution was removed.

**Checkpoint:** retain the conditional-history interpretation; settle the
network stress/loss matching and the physically matched charged-exit state
before combining the modules into a model-to-observation chain. The action
still does not independently fix its physical scale or primordial charge
yield. No source of net X charge is added by this review.

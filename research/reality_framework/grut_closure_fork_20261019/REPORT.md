# GRUT after Candidate 1: deciding the closure fork

19 October 2026. Parent: `8fdd674d5e487dc050d685890e7c6defc35e18b1`.

The review is right about the research priority: Candidate 1's present
benchmark is a conventional, strongly classically bound scalar object in a
weakly coupled EFT. Its calculated loop coefficient provides a model-specific
consistency result; it does not supply a new physical principle. The proposed
closure replacement, however, relies on thermality and curvature claims that
need correction before they can support Option A.

The work completed here confirms the cheap classical channel comparisons,
rechecks the actual slow-roll statement, and gives explicit counterexamples
to the proposed generic curvature falsifier and temperature transfer. The
result supports preserving the existing model and its claim boundaries,
while defining a much more precise target for any cosmological successor.

## Candidate 1 and the value of another loop calculation

Graham's [2001 Q-ball paper](https://arxiv.org/pdf/hep-th/0105009) is prior art
for the spectral approach. It was already cited in the energy packet's
[derivation](../grut_candidate_1_energy_20261019/ENERGY_DERIVATION.md) and
[challenge](../grut_candidate_1_energy_20261019/CHALLENGE.md). The earlier report
does not claim invention of the method. The coupled Candidate 1 operator and
matching differ from that paper's model; its numerical correction cannot be
imported as our result.

The observed loop shift `.00768%` is about a thousand times smaller than the
classical free-carrier binding fraction. A weak-coupling expectation does not
provide a numerical bound on every infrared-enhanced fluctuation sum, but this
benchmark was a poor place to expect a binding-sign surprise. Verification
counts measure executed checks; they are not evidence for a new law or the
existence of the added fields in nature.

The inexpensive follow-up is now executed at the same couplings. Six BVP
solves cover N39944, N39943 and N19972 with both domain and tolerance refinement.
At **tree order throughout**, in units gv:

| Channel difference | Result |
|---|---:|
| E0(39944)-E0(39943)-mX(tree) | -0.199999152115 |
| E0(39944)-2E0(19972) | -3067.397336565 |

The refined energies change by less than `1.4e-9 gv`. The exact classical
identity is `dE/dQ=omega`; the unit finite difference is its integral, here
`.800000847885 gv`. The endpoint approximation is entirely adequate for the
emission sign. Both selected channels are energetically closed at tree order.
These results support deprioritizing a large quantum channel calculation at
this parameter point. They are not substituted for uncomputed quantum energy
differences or a lifetime.

## What v11 and v14 actually establish

The stored v11 scope is a conformal scalar's **fixed-mean-time** two-point
function in **spatially flat** FLRW. For a slowly changing Hubble parameter,
it agrees through first order with the de Sitter kernel at

\[
T_{\rm eff}=\frac{H+\dot H/H}{2\pi}
=\frac{H(1-\epsilon)}{2\pi}.
\]

Our exact power-law countercheck reproduces that restricted result. An
endpoint-anchored kernel used by an instantaneous detector has a nonzero
first-order correction at the same nominal temperature; translating both
proper times also changes the kernel at first order. The stored v14 statement
adds a massless minimally coupled sector only with its specified nonzero-
frequency/Wigner and infrared qualifications. Its full mode calculation was
read, not rerun, in this packet.

[Obadia 2008](https://arxiv.org/pdf/0804.2890) already gives a local conformal
coefficient whose slow-roll expansion agrees with this first-order number,
and explicitly separates that coefficient from actual thermality. This does
not establish priority for every finite-separation or minimally coupled
statement. It does rule out presenting the coefficient alone as a new GRUT
prediction. A publishable novelty claim remains unestablished.

There is also a **curvature matching** obligation. In closed global de Sitter,
`a=cosh(H_Lambda t)/H_Lambda`, the Euclidean state's geodesic temperature is the
standard `H_Lambda/(2 pi)`. Yet `a_double_dot/a_dot=H_Lambda coth(H_Lambda t)`.
Their discrepancy is first order in the spatial curvature fraction at late
times. The flat v11 formula cannot be carried onto the closed saddle without
checking that limit. [DERIVATION.md](DERIVATION.md) gives the explicit test.

The known no-boundary-to-Euclidean-vacuum connection is a conditional matter
state construction, not a theorem that every closed universe is thermal.
[Chen, Lin and Yeom 2018](https://link.springer.com/article/10.1140/epjc/s10052-018-6426-4)
also emphasize the distinction between the familiar short-wavelength result
and curvature-sensitive long wavelengths. Closure alone fixes neither the
observer's restriction nor its thermal spectrum. The required action, state
prescription and observation must still be supplied.

## Curvature sign is not a numerical falsifier

Compactness does not imply positive local spatial curvature: a periodic flat
three-torus is compact and has zero spatial Ricci scalar. Particular
no-boundary saddle choices or global equations may impose stronger conditions.

Even if a chosen cosmology requires k=+1, the proposed falsifier does not
follow. By definition `Omega_K=-1/(a^2 H^2)`. In the closed de Sitter limiting
example,

\[
\Omega_K=-\frac1{e^{2N}-1}.
\]

At N5 it is already `-4.54e-5`; at N10 it is `-2.06e-9`. The spatial
sections remain spheres. For a realistic history, the exact dilution identity
`|Omega_K,f/Omega_K,i|=exp(-2N)(H_i/H_f)^2` adds a reheating/history factor.
Without a restriction or probability distribution for N, closure supplies no
minimum present curvature magnitude. A measurement consistent with flatness
at `1e-3` therefore cannot exclude the entire negative-Omega_K sign class.

[Maldacena 2024](https://arxiv.org/html/2403.10510v2) identifies the harder
problem: the usual semiclassical no-boundary weight heavily favors too little
inflation and excessive curvature. That is a problem for a specific
probability prescription, not a universal lower bound from compactness.
Changing its conditioning or measure requires an additional justified choice.

The formulations must also be distinguished. The
[original global sequestering model](https://arxiv.org/pdf/1406.0711)
requires finite worldvolume and a future collapse under its assumptions.
The [local four-form formulation](https://arxiv.org/pdf/1505.01492) explicitly
permits infinite spacetime volume and a residual integration constant.
Assigning the original global restrictions to every sequestering model is
incorrect. Sequestering changes the global gravitational constraints; it is
not simply a state added to otherwise unchanged Einstein dynamics.

## Observation is already more specific than the proposed bet

The historical Planck-only curvature preference is not a generic present
measurement. For example, the DESI collaboration's
[DR2 table 5](https://arxiv.org/html/2503.14738v3) gives
`Omega_K=.0023 +/-.0011` for DESI+CMB in curved LambdaCDM, and
`.0011 +/-.0013` for DESI+CMB+Pantheon+ when evolving dark energy is also
allowed (68% intervals). Dataset and model assumptions matter. These fits
do not establish a positively curved GRUT cosmology.

The specific CMB-S4 forecast also needs updating: the official
[Fermilab closeout seminar](https://indico.fnal.gov/event/72356/) records its
cancellation in 2025. Forecast precision cannot be treated as a guaranteed
future observing schedule. No new empirical likelihood or Euclid forecast
is calculated in this packet.

## Decision on A, B and C

| Option | What is defensible now | What it still requires |
|---|---|---|
| A: closure as a completed reading | An exposition of explicitly adopted GR/QFT and a named state prescription | The stated new thermality theorem and universal curvature falsifier do not presently follow; empirical content needs a complete model and distribution |
| B: research on closure | A real research direction with sharply stated mathematical choices | A contour/measure or state-selection rule, normalizable perturbations and a calculable observable distinction; existing disputes are not resolved by renaming the premise |
| C: preserve the map | The current result: specified models, conditional derivations and documented failures | Preserve each theorem's actual hypotheses; failure of tested mechanisms is not a theorem excluding every possible responsive-vacuum law |

The recommended current status is C, with Candidate 1 preserved as a numerical
EFT benchmark. Work seeking additional physical content should pursue one
specific question within B: **does an explicit closure condition select a
cosmological quantum state or integration cycle that the supplied dynamics
otherwise leave undetermined?** No such selection rule has been established
in the material inspected here. The explicit research contract is in
[NEXT_QUESTION.md](NEXT_QUESTION.md).

This does not declare the program permanently impossible. It identifies the
missing mathematical choice that another construction would have to resolve.
A state principle can add scientific content by predicting preparations or
cosmic initial conditions; it need not derive new equations to be useful.
It must actually select the state, however. A name or a compact topology does
not make that selection.

## Execution and scope

The finite calculations are in `probe.py` and `results.json`. A missing
mpmath dependency was handled with 60-digit built-in Decimal arithmetic,
without installing a package. The first finite-epsilon grid failed one
asymptotic-ratio criterion at its coarsest point; `results_initial.json`
preserves it. Two additional, smaller epsilon values bring the ratios toward
four and reproduce the predicted second-order coefficient. The refinement
and geometric transfer checks were recorded before their execution.

No no-boundary path integral, cosmological fit, new Q-ball loop spectrum or
quantum lifetime was computed. No independent external validation or
publication novelty certification has occurred. Prior packets remain intact.
Sources and actual read levels are in [SOURCES.md](SOURCES.md).

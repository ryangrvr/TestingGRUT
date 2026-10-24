# Fast conversion products: a source history and a timing-dependent streaming test

23 October 2026 · Codex / AInstein · parent bb211d1

**The review's central concern is correct: hidden conversion creates fast X,
and its production time matters.** At the unchanged reference point, the
worked early 20-MeV source gives a mean post-escape horizon of about 195
comoving pc, whereas a source at z=10^6 gives 8.38 Mpc. A source at 1 keV
already gives 2.37 Mpc. These are conditional transport results, not a
Lyman-alpha exclusion or a completed formation history.

## What was actually calculated

Nine forward/inverse-conversion evolutions cover three grids, time steps
and source-bin widths for each imposed epoch. The rate table is the
parent's fixed Q=200000 ground-object calculation. The same point,
instantaneous scalar-at-rest daughter source and original supplied Y_Q
are retained. Neither the charge yield nor v is retuned. Object profiles
do not evolve, and this calculation is not the unresolved core-formation
source or a quantum hot-core lifetime.

An independent optical-depth integral gives the first-conversion CDF.
The coupled population includes inverse conversion. Positive particles
which redshift below its kinematic threshold enter a receiver carrying
their exact comoving momentum and escape time. Collisionless horizons
are then calculated independently of the reactive energy grid.

## The finite-Q energy correction

At this object, **omega/m_X=0.861195131**, whereas
**E/(Qm_X)=0.923316135**. Binding per particle does not by itself
give the derivative dE/dQ. The cold outgoing kinetic energy is
K_out/m_X=2Omega=1.722390262; gamma_out=2.72239
and p_out/m_X=2.532076. For the imposed warm
daughter input, K_out/m_X=1.947135133
and p_out/m_X=2.772292. The outgoing
particle is relativistic at this point. A nonzero kinematic allowance
alone is not a statement about its amplitude in every possible model.

## First conversions occur promptly, then inverse reactions matter

The initial rates Gamma_conversion/H are approximately 3.23e6 at the
20-MeV epoch, 115 at 1 keV, and 27.0 at z=10^6. For the last case,
half the initial anti-X have first-converted by **z approximately
974487**. The median exit
from the reactive band occurs later, near z approximately
890452.
Those are different source times. Forward events, inverse events and
net anti-X abundance are recorded separately; repeated conversions are
not counted as newly captured particles.

At 20 MeV first conversions occur well inside one ordinary transport
step, so their timing is obtained from the independent survival integral.
Its 32-to-64-point quadrature check differs by less than 1.3e-15 in the
recorded probes. This checks the supplied rate interpolation, not the
unresolved accuracy of every narrow scattering resonance.

## The resulting scales

All distances below are comoving and number weighted within the converted
population. The upper bound adds the largest possible traveled prefix
from source to escape, using the speed of light. It bounds traveled path,
not a matter-power cutoff or a lower bound on net spatial displacement.

| Imposed source epoch | Mean horizon to z=5 | Mean horizon to today | Total traveled-path mean upper bound to today |
|---|---:|---:|---:|
| 20 MeV | 194.46 pc | 194.95 pc | 195.34 pc |
| 1 keV | 2.35546 Mpc | 2.36919 Mpc | 2.38105 Mpc |
| z = 10^6 | 8.31841 Mpc | 8.37832 Mpc | 8.44154 Mpc |

Across the three refinements the largest relative horizon spread is
0.0653% (the early case). The pre-escape causal allowance is less than
one percent of the mean horizon in these examples. At the 20-MeV source,
the converted horizon is small compared with Mpc scales. The late sources
require a genuine mixed-matter perturbation and structure-data test.
Elastic deflection, angular conversion kernels and source correlations
prevent identifying these path lengths directly with a density cutoff.

## Energy fractions are not warm rest-mass fractions

The imposed source initially has anti-X rest mass at **3.1751%** of the
old calibrated matter budget. Complete net conversion releases
**5.4688%**
of that budget from core chemical energy, approximately
996.512 GeV per reference object.
The fast energy is a real issue; it is not eliminated by calling it hidden.

After that kinetic energy redshifts away, the converted X have
**3.3588%**
of remaining rest matter; **all free X have
6.7176%**.
The initial positive daughters are a second free component with their own
smaller horizon. Which part is operationally warm depends on the source
epoch and the resulting perturbation spectrum. The roughly 7% birth-energy
estimate is not a 7% warm rest-mass fraction.

This completion also leaves only
**94.5312% of the old
late rest-matter calibration**, within the frozen-profile approximation.
The old abundance and detector example therefore cannot simply be carried
over. Core charge transfer is about 6.26%; evolving profiles, the actual
formation source and updated abundance must be calculated before a new
population forecast is certified.

## A numerical limitation was found and excluded

An attempted active-grid estimate of total path underestimated an independent
no-reaction control by **8.10%**. Its recorded fields are retained as rejected
diagnostics. They are not used for the horizons or physical conclusions above.
The comoving-momentum receiver, independent analytic free-flight integral
and causal prefix bound provide the released result instead.

The **8/8 accepted validation groups pass**: independent
radiation-era quadrature, extremely cold momentum arithmetic, number/charge/
energy ledgers, three refinements, input identity, first-event quadrature,
mass accounting and explicit fencing of the rejected path estimator.
These are author checks, not independent replication or empirical validation.
The rejected diagnostic remains visibly failed rather than absorbed into
the pass count.

The actual soft formation source, general hot-core leakage and vacuum density
remain unresolved. The next observational calculation must evolve mixed-matter
perturbations with this specified source and the required angular/scattering
physics. No thermal-WDM mass bound is relabeled as an exclusion here.

![First conversion timing and the conditional converted-product horizons](conversion_streaming.png)

[Equations](MATHEMATICS.md) · [results](STREAMING_RESULTS.json) ·
[first-event history](FIRST_CONVERSION_RESULTS.json) ·
[verification](VERIFICATION_RESULTS.json) ·
[rejected diagnostic](REJECTED_PATH_DIAGNOSTIC.json) ·
[Claude handoff](CLAUDE_HANDOFF.md) · [reproduction](REPRODUCE.md) ·
[sources](SOURCES.md).

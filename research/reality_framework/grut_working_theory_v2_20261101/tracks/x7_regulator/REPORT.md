# X7 · Spatial cutoff and coherent H55 overlap

6 November 2026 · TestingGRUT · [Model and derivations](MODEL.md)

The signed auxiliary response controls the sharp spatial edge in
the one-axis expanding slab. With fixed auxiliary masses, the
refined slab meets the 1% numerical criterion for accumulated and
instantaneous work. Changing the auxiliary masses still changes the
instantaneous work by 11.5% at the wider cutoff. The H55 homogeneous
restart also fails the tested transport windows. The source and
transport gates for an H55 history remain open.

## Spatial source

The unsubtracted missing-neighbor term grew approximately as $Q$.
The signed term fell as $Q^{-3.00}$ to $Q^{-3.03}$ across the three
harmonics and two tested scale factors. The signed static response
at the largest cutoff with fine transverse bins agrees with the
independent finite-periodic result within 0.28%. A coarser radial rule
at that cutoff leaves errors up to 5.72%. The separate occupied-to-empty
spectral evaluation agrees with the finite frequency-dependent kernel
within $1.6\times10^{-16}$ absolute in the tested units.
[Static and frequency response](results/static.json)

The matched energy ledger compares integrated reciprocal work to
stress change plus expansion work. The signed homogeneous slab
recovers the independent helicity ODE to $2.80\times10^{-7}$ maximum
absolute occupation error at the finest step. Its signed energy and
pressure agree within $3.33\times10^{-6}$ relative. The volume and
mass variations, static vacuum, massless decoupling and fixed-scale
reversal also pass. Across the expanding runs the largest total
pressure-work budget residual is $9.21\times10^{-8}$ and the
largest comoving momentum residual is $2.60\times10^{-11}$.
[Controls](results/controls.json)

| Matched-work comparison | Peak-scaled instantaneous error | Cumulative history error | Terminal transfer error |
|---|---:|---:|---:|
| Initial wider cylinder, (12,16) to (16,22) | 0.0865% | 0.0124% | 0.00665% |
| Initial transverse midpoint width 0.5 to 0.25 | **1.246%** | **1.052%** | **1.039%** |
| Finer transverse midpoint width 0.25 to 0.125 | 0.297% | 0.259% | 0.256% |
| Fine midpoint to independent panel Gaussian rule | 0.0985% | 0.0857% | 0.0850% |
| Gaussian rule 8 to 12 nodes per panel at (16,22) | 0.0928% | 0.00277% | 0.00277% |
| Wider cylinder with the 12-node Gaussian rule | 0.0829% | 0.0126% | 0.00664% |
| Timestep halved in the original cylinder | 0.00312% | 0.000219% | 0.0000431% |
| Heavier auxiliary pair, original cylinder | **11.61%** | **1.805%** | 0.996% |
| Heavier auxiliary pair, wider Gaussian cylinder | **11.53%** | **1.795%** | 0.990% |

Peak scaling uses the reference curve's maximum absolute instantaneous
work. Cumulative scaling uses the maximum absolute cumulative transfer.
The independent longitudinal cutoff, transverse cutoff, scalar
collocation and scalar bandwidth comparisons also meet 1%.
[Original comparisons](results/evolution.json) ·
[Radial follow-up](results/radial.json)

The compact pulse leaves a weighted auxiliary excitation energy
of $9.71\times10^{-5}$ of the physical excitation energy at the
finer timestep. This control tests real excitations after the drive
ends. The matched vacuum polarization is separately mass dependent.
For harmonic one at $a=1$ and $\omega=2.6$, the difference between
the two auxiliary pairs after removing the cutoff is
$-1.097\times10^{-3}$ in the quadratic response kernel. At the
largest finite cylinder the edge difference is $1.22\times10^{-7}$
and the total kernel difference agrees with that unbounded value
within $1.1\times10^{-4}$ relative. The kernel calculation isolates
a finite virtual effect. The live work comparison also includes
nonlinear response, expansion and preparation.
[Compact pulse](results/pulse.json) · [Frequency response](results/static.json)

The auxiliary pairs $(3.2,4.8)$ and $(4.8,7.2)$ and the control's
scalar profile are choices for the synthetic slab. The result
establishes discretization control at fixed auxiliary masses. The
regulator-independent source still needs higher-derivative matching
or a controlled auxiliary-mass limit with the cutoff increased
independently. Both need a common action for force and stress. A
spatial H55 run also needs its physical initial profiles.

## H55 interface

The restart reads the archived 5,760-bin, four-species, two-helicity
Bloch state at $t=2.9980937H_i^{-1}$. It evolves the existing scalar
and expansion equations for another $0.12H_i^{-1}$. The second
timestep is half the first. A separate sensitivity applies X5's
instantaneous-basis vacuum decay at the inherited widths. It uses
the same saved input and deposits the lost occupation energy in the
visible radiation ledger. Spectator damping stays in both runs.

The archived state has a dressed energy gap of $9.76\times10^{-4}$,
a free projection pressure gap of $-7.57\times10^{-5}$ of the
total density, and normalized scalar work of $0.10835H_i$.
The restart reproduces all three. For testing later times, work is
normalized by $\rho H_i$ and the energy and pressure gaps by $\rho$.
Each of the four runs has 916 sampled times and **zero** simultaneous
passes of the $10^{-6}$ gate. Every tested $0.02H_i^{-1}$ window has
a violation. On the finer runs the smallest window-maximum gap is
0.890 with coherent evolution and 0.883 with vacuum damping.
The largest time-step difference in the work history is 0.00550.
The total pressure-work residual stays below $4.46\times10^{-12}$.
[Coherent overlap](results/overlap.json)

Reallocating the known scalar kinetic subtraction to the fermion
stress gives initial matched work of $-0.08380H_i$ per unit total
density. The corresponding free projection also fails at every
sample, with energy and pressure compared at fixed scalar velocity.
The vacuum-damping sensitivity puts 0.137% of total energy into
visible radiation by the end of the interval. The archived 50.96%
is the instantaneous physical-mass fermion energy fraction at the
initial slice.

In the exact homogeneous moment equations, one coherence component
carries the pressure correction and the other carries the work.
Reflecting the work component of the archived state keeps its
occupations and instantaneous stress fixed while reversing the
work. The two admissible instantaneous covariances differ in
fermion energy derivative by $0.21670\,\rho H_i$. The archive
reconstruction and energy continuity identity agree at roundoff.
This establishes the missing state information for a transport
overlap using this Hamiltonian.

The separate fixed-$a$, drive-off control checks a sufficient
pressure guarantee for arbitrary phases under vacuum decay.
At half the initial parent energy its pressure envelope is still
$6.44\times10^{-3}$ of the initial total density. That half-energy
time is $44.35H_i^{-1}$, compared with the initial Hubble time
$6.06H_i^{-1}$. This conditional envelope measures robustness to
phase choice and gives no live H55 timing prediction.

The earlier 1.28% terminal N4 count refinement miss and the
inherited strict result of zero accepted points out of 135 are
recorded in the evidence. X7 tests standard semiclassical dynamics
in the matter candidate. A GRUT-specific prediction requires a
derived selection or coupling law with a discriminating observable.
[Verification](results/verification.json) · [Replay](REPLAY.md)

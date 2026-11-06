# X4 · Repairing the H55 fermion interface

5 November 2026 · GRUT M2 · $H_i=5.5\times10^{10}$ GeV

The original H55 source cannot supply a physical free-Majorana distribution
with both of its recorded N1 moments. This survives correcting the time
labels and allowing arbitrary injection shapes with the original source
moments, free streaming and vacuum decay. A new kinetic implementation
closes the full energy ledger with bounded occupations and spectator decay
for two explicit preparations. A separate quantum calculation on the actual
H55 mean waveform identifies a substantial N1 production discrepancy.
Scientific acceptance remains false.

## Source of the mismatch

The passive 64³ scalar replay produced 9,204 timestep records. The archived
scale factor, Hubble rate and total energy density agree within
$5.0\times10^{-10}$ relative at the saved times. Rebuilding the energy and
number recurrence from each trace row agrees to $3.7\times10^{-16}$ of the
peak. Both moments are physical densities divided by the same $v^2$ in the
scalar engine. The kinetic integrals sum two helicities with
$E(p)=\sqrt{p^2+M^2}$ and include rest energy.

The archived row prints energy before the fluid update and count after it.
At the last archived time, $t=2.99809375H_i^{-1}$, aligning the two changes
the number by about 0.0114%. The source evolves energy with $a^{-4}$ and
number with $a^{-3}$ between injections, and assigns number at injection
from a prescribed energy per particle. Its separate blocking energy does
not impose a cap on the resulting momentum distribution.

| N1 calculation at the archived interface | $\rho_N/\rho_{\min}(n_N)$ |
|---|---:|
| Aligned native moments | 0.43642 |
| Massive free streaming of each cohort at its native mean birth energy | 0.44457 |
| Same cohort control with the perturbative floor removed | 0.53692 |
| Same energy requests with stationary two-body energies on the floor branch | 0.23942 |

The cohort controls use the actual held injection ledger. They change
propagation or injection kinematics without feeding the changes back into
the scalar solution. The stationary two-body control increases the
required count because the floor often used the much larger
amplitude-dependent energy per particle. In that control N4 also falls
below the Pauli minimum.

### Bound without a spectral-shape assumption

For each injection at scale factor $a_b$, define $r=a_b/a_t\leq1$ for the
target scale factor $a_t$. Convexity gives
$E(rp)\leq rE(p)+(1-r)M$. Time dilation bounds the survival probability
below by $\exp[-\Gamma_{\rm rest}(t-t_b)]$ for vacuum decay. If
$\Delta n_b,\Delta\rho_b$ denote the physical birth moments, all free
streaming spectra with those moments obey

$$
n_t\geq n_{\rm lo}
=\sum_b r^3 \Delta n_b e^{-\Gamma_{\rm rest}(t-t_b)},\qquad
\rho_t\leq \rho_{\rm hi}
=\sum_b r^3 \big[r\Delta\rho_b+(1-r)M\Delta n_b\big].
$$

The energy bound allows zero decay, while the number bound allows the
fastest vacuum decay at every momentum. This makes the combination
generous. Because the filled-sphere minimum is increasing in number,
the complete N1 injection ledger gives
$\rho_t/\rho_{\min}(n_t)\leq \rho_{\rm hi}/\rho_{\min}(n_{\rm lo})
=0.83074$. This excludes a repair using only a different injection shape,
free streaming and vacuum decay at the fixed birth moments. The
floor-removed version of this particular bound is inconclusive at 1.09942.
Its result is distinct from the floor-removed cohort diagnostic above.

There is a separate kinematic problem with interpreting the old unblocked
floor as perturbative decays of parents at rest. Under the inherited generous
stationary-parent support, $E_{\max}=(m_{\rm parent}^2+M^2)/(2m_{\rm parent})$
even with massless visible daughters. Fill every mode in the final comoving
support ball, permit occupation one in both helicities, and replace every
vacuum loss at the maximal rest rate. With
$p_{\max}=\sqrt{E_{\max}^2-M^2}$ and
$n_{\rm cap}=p_{\max}^3/(3\pi^2)$, an upper bound on injected comoving energy is

$$
E_{\max}n_{\rm cap}\left[a_t^3+\Gamma_{\rm rest}\int_{t_0}^{t} a^3\,dt'\right].
$$

The actual floor requests exceed this generous bound by factors $9.92\times
10^3$ for N1 and $2.99\times10^4$ for N4. This diagnoses the attempted
stationary-parent interpretation. It leaves coherent high-amplitude
production and moving-parent kinematics as distinct mechanisms.

## Spectral evolution with the spectator included

The implementation computes number, energy and pressure from the same
occupation grid. At fixed scale factor the scalar-to-singlet and
singlet-to-visible exchanges conserve energy within $3.2\times10^{-15}$
relative across the eight cases. The maximum pressure-work error is
$4.04\times10^{-4}$ in the kinetic histories.

The spectator obeys
$\dot\rho_s=-(3H+\Gamma_s)\rho_s$ and
$\dot\rho_{\rm DR}=-4H\rho_{\rm DR}+\Gamma_s\rho_s$.
The evolution treats it as pressureless and uses the inherited scale rule
$\Gamma_s/H_i=5.5\times10^{-7}$. Six independent constant-$H$ ODE
controls agree within $1.7\times10^{-14}$ relative. The visible bath
and the dark daughters stay separate in the Friedmann ledger.

| Preparation | Late total $K$ | Late visible $K$ | Full initial energy gap |
|---|---:|---:|---:|
| Actual H55 end slice with energy-matched half-filled spectra | 11.03207 | 11.02318 | $2.1\times10^{-15}$ |
| Same case with the spectator omitted | 11.02309 | 11.02309 | $3.24\times10^{-5}$ |
| Reduced harmonic evolution started at the first waterfall maximum | 12.86983 | 12.85827 | $8.9\times10^{-16}$ |

Here $K=a^2H$ and visible $K$ uses the visible energy alone. Momentum
doubling and timestep halving change both late values by at most 0.094%
in the matched case and 0.204% in the harmonic case. The latter can
restart from its full occupations at the archived interface with an
energy gap of $7.8\times10^{-14}$. Its late $K$ differs by 0.0082% from
the uninterrupted run.

The matched spectra lower the initialized N1 count from
$6.9232\times10^9$ to $3.1256\times10^9 H_i^3$. The harmonic preparation
returns the 35.50% of source energy already in fermions at
$t=0.24840625H_i^{-1}$ to their respective parents, sets occupations to
zero, and replaces the later lattice dynamics by cycle-averaged
condensates. Its restart preserves the occupations of this reduced
model. The actual scalar source still has its 1.33043% budget error.
Both cases keep the common 0.05 eV seesaw input, the inherited
hard/thermal split prescription and the specified spectator approximation.

The late dark-to-visible energy ratios are 0.001612 and 0.001800. With
standard visible entropy evolution from $g_*=g_{*s}=106.75$, these would
give $\Delta N_{\rm eff}=0.00461$ and 0.00514, respectively. Those values
belong to the two stated preparations. [Entropy conversion](https://doi.org/10.1103/PhysRevD.108.063030)

## Quantum evolution on the actual mean waveform

The collisionless Majorana helicity equations use the actual recorded
spatial means of $\sigma,\chi$ and their velocities. Each computed
occupation is bounded by one. The static control has occupation below
$1.6\times10^{-32}$. The refined time-dependent control agrees with an
independent unrotated-basis ODE at $1.5\times10^{-10}$ in occupation.
Five H55 grids control momentum spacing, cutoff and integration step.

| Species | Fine-grid dressed energy / old fluid energy at the interface |
|---|---:|
| N1 | 0.51831 |
| N2 | 1.05379 |
| N3 | 0.98473 |
| N4 | 0.94881 |

Across the five controls the N1 ratio is 0.5130–0.5204. At the terminal
slice, doubling momentum resolution from 1,440 to 2,880 changes N1 by
0.38% and N4 by 1.27%. Halving the integration step and widening the
cutoff at fixed spacing change the N1 energy by less than 0.002%.
The two vacuum projections give close terminal moments on the tested
waveform. These are free-particle diagnostics of a prescribed,
collisionless, spatially averaged background. Inhomogeneous modes,
collisions, renormalized stress and feedback can change the result.

## Decision

The source replay and bound close the old moment-matching route under its
own particle assumptions. The numerical ledger and exact restart supply
a viable interface for a new explicit spectral state. The physical H55
history still needs a self-consistent production source. The next
calculation must obtain particle production and scalar work from the
same state on a resolved background. At a finite regulator the
collisionless identity is

$$
\frac{d}{dt}{\rm Tr}(H_F C)={\rm Tr}(\dot H_F C)
\quad\text{when}\quad i\dot C=[H_F,C].
$$

The scalar force must carry the opposite contribution from the
scalar-dependent part of $\dot H_F$. A renormalized implementation
needs a common subtraction and counterterm prescription for the scalar
force and stress, plus a controlled accounting of expansion and decays.
Only that coupled history can decide if an actual H55 slice passes the
source budget and carries a physical spectrum into the late solver.

The strict scan remains zero accepted keys of 135. There is no new
release result, object population, cosmological root or validated GRUT
observable. All 39 evidence checks are internal. The archived code and
results are available for independent replay.

[Methods](METHODS.md) · [Replay](REPLAY.md) · [Handoff](HANDOFF.md) · [Evidence checks](results/verification.json)

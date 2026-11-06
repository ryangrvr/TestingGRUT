# X3: the actual H55 continuation and the fermion-moment obstruction

5 November 2026 · resumed from `9f71d34f8a768e4ec7bc60a46a65a833a818c1de`.

**Result:** the actually evolved lower-scale scalar history gives a late
reheating-network \(K=a^2H\) of **11.0312**, compared with **7.01656** for the
old H55 proxy under the same new spectral prescription: **+57.22%**. The
conditional network result is numerically stable in the controls below.
However, the native N1 energy and number density cannot both be represented
by a free, two-helicity Fermi distribution. A faithful native-moment splice
therefore **fails**. No cosmological point or relic population is accepted.

This finishes the registered X3 investigation as a conditional numerical result
and a closure failure, not as a successful matter model. It does not exclude
M2's field content, all possible initial states, or GRUT as a framework.

## 1. What was recovered and what actually ran

The previous chat's public share page could not be loaded. The current state was
recovered from the repository and recent conversation context. The authoritative
sequence is:

| Commit | Meaning |
|---|---|
| `8ab5829a` | Completed strict M2-X scan; fresh H55 scalar control |
| `98ccd20d` | Emergence-and-interaction synthesis; no new cosmological run |
| `9f71d34f` | Registered actual-history continuation; no kinetic result yet |
| `d845c6d` | X3 execution protocol, including the newly noticed moment issue, before new kinetic runs |

Read the original [scope](SCOPE.md), [execution protocol](RUN_PROTOCOL.md),
[replay commands](REPLAY.md), [machine-readable summary](results/summary.json)
and [verification](results/verification.json). The [manifest](MANIFEST.json)
pins package bytes and inherited source identities. All new calculations and checks
are same-author work. No outside review or independent-agent review is claimed.

Executed here: four transition-table replays; seven full dephased kinetic
continuations; exact and numerical Fermi-moment feasibility checks; nine
homogeneous line controls; nine full-3D **single-mode** controls; population
conservation identities. Not executed: a new scalar lattice, a full H55 Dirac
sea with backreaction, a recalibrated H55 drain, a release/root scan on the new
history, or a cosmological nugget population.

## 2. A necessary condition the native moments fail

The archived raw diagnostic stores `rN1..4` and `nN1..4`. The source's units are
energy density divided by \(v^2\), and number density divided by \(v^2\),
respectively. Multiplying by the recorded \(v^2\) gives the Hubble-unit inputs
of the kinetic engine. Its `species_energy` routine sums **two helicities**
with \(E(p)=\sqrt{p^2+M^2}\); this is the convention tested here.

For \(0\le f_h(\mathbf p)\le f_{\max}\), the minimum energy at fixed number
is obtained by filling the lowest-energy states first. With \(g=2\),

\[
n=\frac{g f_{\max}k_F^3}{6\pi^2},\qquad
\rho\ge\rho_{\min}(n,M)=\frac{g f_{\max}}{2\pi^2}
\int_0^{k_F}p^2\sqrt{p^2+M^2}\,dp.
\]

The integral is

\[
\frac18\left[k_F E_F(2k_F^2+M^2)
-M^4\operatorname{asinh}(k_F/M)\right].
\]

Proof: moving any occupation from a higher-energy occupied state to an available
lower-energy state preserves number and reduces energy. Repeating the exchange
fills a sphere. Anisotropic occupations cannot lower this minimum. Within the
same local free-gas approximation, spatial concentration cannot lower it either:

\[
\frac{d\rho_{\min}}{dn}=E_F,\qquad
\frac{d^2\rho_{\min}}{dn^2}=\frac{k_F^2}{3nE_F}>0,
\]

so Jensen's inequality gives the uniform-density lower bound for a spatial
average. This is a standard Fermi-gas variational constraint, not a new GRUT law
or a theorem about a general interacting stress tensor.

At the final recorded H55 time, \(t=2.9980937\), the measured ratios are:

| Species | \(\rho/\rho_{\min}\), physical Pauli cap 1 | \(\rho/\rho_{\min}\), dephased cap 1/2 |
|---|---:|---:|
| N1 | **0.436485** | **0.346667** |
| N2 | 0.993325 | 0.789531 |
| N3 | 1.004889 | 0.797939 |
| N4 | 1.021774 | 0.811142 |

N1 has \(\rho_1=1.3385157\times10^{13}\) and
\(n_1=6.9223981\times10^9\), in \(H_i^4\) and \(H_i^3\) units. Its physical
minimum energy at that number is \(3.0665822\times10^{13}\). At fixed energy,
the number must fall by **at least 46.34%** to satisfy the physical cap, or
**54.85%** under the additional half-occupation prescription.

The 0.67% N2 physical-bound deficit is small compared with the source's maximum
scalar budget diagnostic, **1.33043%**. It is not a robust independent physical
exclusion. The much larger N1 deficit is the decisive obstruction. The printed
fractions and fluid densities also differ at a few parts per million because
of output precision and the source's update order; this cannot repair N1.

Checks: direct quadrature agrees with the closed integral to floating-point
precision. A separate 4000-bin linear program reproduces the minimum to
\(1.45\times10^{-7}\) relative. Seven saved times were tested. See
[moment_feasibility.json](results/moment_feasibility.json).

N number is not a protected charge in this model. A physical reaction could
change it. But no finite-time reaction or alternative particle dictionary has
been calculated to connect these incompatible moments to the supplied spectrum.
Replacing the number at the interface must be labelled a new preparation rule.

## 3. The conditional kinetic continuation

The baseline uses the actual H55 scalar energy at its final **recorded** time,
not an extrapolation to \(t=3\). It keeps the old kinetic reference's species
energy shares and SM hard/thermal split to isolate the history change. A second
run uses H55's native species energy shares. All start from an explicitly
supplied spectrum: occupation 1/2 in both helicities below an edge, with one
partially occupied bin so each energy is matched exactly. The dilute variant
uses occupation 1/4. The resulting number is recorded, never silently clipped.

The unchanged dephased production/collision/expansion engine then runs to its
original stopping condition. The original four transition tables replay
**exactly**. The timestep control changes only three adaptive step ceilings in
an isolated source copy. In these runs the seesaw rate input is the inherited
common \(\widetilde m_j=0.05\) eV for all four species; it does **not** silently
substitute the later card's \(0.09122\) eV N2 value. That card-to-network
matching remains a separate model dependency.

| Run | Momentum bins | Late K | Change from actual baseline |
|---|---:|---:|---:|
| Old H55 proxy, matched spectral rule | 6400 | 7.01655683 | −36.39% |
| Actual H55, coarse grid | 3200 | 11.02834213 | −0.02623% |
| **Actual H55, baseline** | **6400** | **11.03123521** | — |
| Actual H55, refined grid | 12800 | 11.03125904 | +0.000216% |
| Actual H55, halved timestep ceilings | 6400 | 11.01701376 | −0.12892% |
| Actual H55, native species energy shares | 6400 | 11.03113764 | −0.000884% |
| Actual H55, dilute supplied spectrum | 6400 | 10.81040152 | −2.00189% |

The matched proxy is within 0.004% of its archived K=7.01683166. The baseline
ends at \(t=252084.26\), with zero condensate fraction at printed precision and
the original N-energy stopping condition satisfied. Its last-decade K drift is
\(3.59\times10^{-5}\) relative. The reported SM-half time is \(10164.20\), with
\(T_{\rm RH}=1.30763\times10^{12}\) GeV. Threshold times are sampled at output
cadence; the K comparison is the better-resolved diagnostic.

**Why this matters:** in a radiation-dominated network with fixed relativistic
degrees of freedom, \(a^3s\propto K^{3/2}\). Holding the produced comoving charge
fixed, the matched K shift alone would change its entropy-normalized yield by

\[
\frac{Y_{Q,\rm actual}}{Y_{Q,\rm proxy}}
=\left(\frac{11.031235}{7.016557}\right)^{-3/2}=0.50728.
\]

This is a conditional scaling consequence, **not** a recalculated relic yield:
the actual history also changes the release forcing, and that calculation did
not run. There is no license to move or accept an old root using K alone.

## 4. Energy conservation: the passed and failed checks

Within the **represented reheating network**, each supplied species energy is
matched to below \(10^{-15}\) relative, and the engine's initial Friedmann
energy agrees with that network to at worst \(4.44\times10^{-11}\).

The **whole archived scalar source** contains a spectator energy fraction
\(3.23821\times10^{-5}\), which this inherited kinetic engine does not carry.
With the printed-column residual, the actual full-source energy gap is
\(3.33907\times10^{-5}\); the Hubble jump is about \(-1.67\times10^{-5}\).
That exceeds the registered \(10^{-6}\) full-splice tolerance. It is explicitly
recorded as a failure. The original 1.33043% scalar budget error is also retained.
Neither is renormalized away.

A complete energy history must evolve the spectator and its dark-radiation
daughters under their stated dynamics. The present output is a conditional
reheating-network K, not the full late-time universe. A conservation check on
selected components cannot certify an omitted component.

## 5. Controlled quantum tests and their limit

The inherited homogeneous line algorithm was compared with an independently
integrated two-component Schrödinger equation using DOP853 and eigenvectors
obtained by direct matrix diagonalization. At amplitude A=10, halving dt from
0.004 to 0.002 to 0.001 reduces the maximum sampled occupation error from
\(3.40\times10^{-5}\) to \(8.50\times10^{-6}\) to \(2.13\times10^{-6}\).
Halving momentum spacing changes the integrated energy by less than \(10^{-9}\)
relative in this resolved test. A=42.3 was also tested. Maximum norm error across
the line controls is \(2.69\times10^{-12}\).

Static-vacuum split-step production is **reported separately**, not subtracted
to make the result pass. At A=10 its integrated artifact falls from
\(7.32\times10^{-4}\) to \(4.57\times10^{-5}\) to \(2.85\times10^{-6}\) as dt
halves. The independent static solution has zero production.

A new **full 3D single-mode** calculation includes all three momentum components
and a spatially varying mass phase. It evolves the four-component Hamiltonian

\[
H=\begin{pmatrix}\boldsymbol\sigma\!\cdot\!\mathbf k&M e^{i\theta(\mathbf x,t)}I_2\\
M e^{-i\theta(\mathbf x,t)}I_2&-\boldsymbol\sigma\!\cdot\!\mathbf k\end{pmatrix},
\]

using exact kinetic and local mass factors. In a \(2\pi\) periodic box, M=1,
\(\theta=1.2\sin^2(\pi t/2)[1+0.3\cos x+0.4\cos y+0.2\cos z]\),
and one initial negative-energy mode at \((1,0,1)\), the final positive-energy
projection is 0.22170134 at dt=0.01. Its change between 16³, 24³ and 32³ is
below \(6\times10^{-13}\) relative. Halving dt changes it by 0.00291%.
The norm error stays below \(3.17\times10^{-13}\). A homogeneous control agrees
with an independent four-component ODE; a static control's projection artifact
falls from \(9.92\times10^{-11}\) to \(6.20\times10^{-12}\) on halving dt.

This establishes an implemented, resolved quantum-mode propagation check. It
does **not** establish a produced-particle density, a filled quantum vacuum,
vacuum-stress renormalization, backreaction, the H55 drain, or cosmological
continuum convergence. An initial negative-energy mode is not a whole Dirac sea.

The inherited cosmological method remains a 1+1D ray with a 3D radial weight.
That weight is exact in the homogeneous reduction; it does not add transverse
streaming or helicity mixing. Stochastic fermion methods are established in the
literature, but their existence does not certify this dimensional reduction:
[Borsányi–Hindmarsh, arXiv:0809.4711](https://arxiv.org/abs/0809.4711).

For the reported H84 hot-spot momentum \(k_{\max}\simeq1.2\times10^5H_i\) and
box \(L=0.42/H_i\), the direct uniform-grid Nyquist bound is about \(1.6\times10^4\)
cells **per axis**, over \(4.1\times10^{12}\) sites. One complex128 four-spinor
alone needs about 264 TB, before replicas, work arrays or refinement. This is a
resource estimate for that method, not a proof that another controlled method
is impossible. H55 calibration has not been performed.

## 6. Population equations without an invented population

Use separate comoving core distributions \(F_Q(\varepsilon,\mathbf P)\), free
\(f_+,f_-\), a coherent charge/energy ledger, and neutral reservoirs. After the
heavy charged fields are gone, a restricted charge ledger is

\[
Q_{\rm total}=\sum_Q QN_Q+N_+-N_-+Q_{\rm coherent}.
\]

If Y or Sigma survives, its corresponding charge must also be included. A
neutral decay is not antiparticle erasure:

\[
\bar X+C_Q\rightleftarrows X+C_{Q-2},\qquad
s_{\rm neutral}\to X+\bar X.
\]

Both reactions conserve charge. The first transfers the exact core mass
difference into particle kinetic energy, with recoil included in an exact
four-momentum treatment. In the heavy, stationary-core approximation,
\(E_{+,out}=E_{-,in}+M(Q)-M(Q-2)\). Its inverse requires the corresponding
incoming energy. Replacing this difference by \(2\,dM/dQ\) is a large-Q
approximation, not a new conservation law.

The code checks 23 charge-neutral reaction vectors and 11 exact synthetic
event-energy balances. It explicitly exhibits different populations with the
same total charge. These are algebraic tests, not formation simulations.
An actual abundance still requires formation, transport, survival and matched
forward/reverse kernels. A stationary solution and a root in an abundance
equation do not supply that data.

## 7. Closure decision

| Requirement | Outcome |
|---|---|
| Actually evolved lower-scale scalar input | Used, with original H84 drain and error retained |
| Conditional kinetic evolution and K controls | Completed; numerical changes below the registered 5% flag |
| Native energy-and-number-preserving particle splice | **Fails for N1** |
| Full-source energy continuity at 1e-6 | **Fails in the inherited reduced network** |
| Controlled H55 microscopic fermion drain | **Unresolved** |
| Actual core population and new accepted root | **Not established** |

The earlier blind result remains **0/135 strict M2-X passes**. This work adds no
accepted point. It closes the attempted faithful splice of this saved instance;
it preserves an informative conditional network calculation. Reopening it needs
a new, declared closure that supplies physically representable spectra and
tracks every relevant reservoir, followed by a fresh scalar-to-kinetic history.
More roots on the old proxy would not address the failure found here.

All **64 evidence checks** pass. Their meaning is that the archived calculations
and the stated failures are internally consistent; the scientific acceptance
flag remains **false**. Replay and source files are available for external review.

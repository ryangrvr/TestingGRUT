# ALP conversion remains unresolved at the frozen benchmark

The homogeneous and linear benchmarks pass. All eight nonlinear runs
transfer energy out of the rotating mode and miss the late particle
gate. Increasing the spacing resolution from $48^3$ to $64^3$
changes the late quadratic count by **6.72%**. Even at $64^3$,
**4.70%** of that count lies near the spectral cutoff and the
anharmonic energy fraction reaches **5.26%**. The late physical
$C$ and photon coupling band remain open.

The physical hierarchy and preparation introduce a separate gap.
For the fixed idealized radiation history, $m_a/H_*=3.60\times10^5$.
Linear growth reaches the registered backreaction warning at
$a\simeq0.448$ for the prescribed classical fluctuations and
$a\simeq0.523$ for the vacuum preparation. Both precede the
nonlinear runs' $a_i=0.6$. The slowest nonlinear run has
$m_a/H_*=200$. Its yield cannot calibrate the physical point.

[Protocol](PROTOCOL.md) ·
[Conversion decision and spectra](results/conversion_decision.json) ·
[Physical linear histories](results/physical_linear.json) ·
[Retrospective ORGAN check](results/organ_retrospective.json)

## Charge normalization and self-consistency

The frozen EFT uses the constant-mass cosine, a constant $f_a$,
and $c_W=-1$ at $\mu=130$ GeV with the other matching coefficients
zero. After electroweak freeze-out, radiation entropy is conserved,
and the scalar is a spectator with no further SM charge exchange
in this model. The mass torque is retained. The calculation uses
$q_\theta=4/\pi$, the early dimensionless rotational action per
$2\pi$ in units of $f_a^2m_a$. It fixes $a=1$ at the adiabatic
separatrix, where the rotational action divided by $2\pi$ is
$4f_a^2m_a a^3/\pi$.
The corresponding libration action is twice as large. This gives
$C=2$ in the slow homogeneous trapping limit. The full ODE follows
the phase through the separatrix to measure finite-rate deviations.
[Homogeneous kinetic misalignment](https://arxiv.org/abs/1910.14152) ·
[Action derivation and fragmentation](https://arxiv.org/abs/2206.14259)

The baryon abundance fixes $Y_\theta/f_a^2=Y_B/(c_BT_{\rm EW}^2)$,
so $s_*=q_\theta m_a c_BT_{\rm EW}^2/Y_B$. Thus the time and scale
factor at which the fixed early angular velocity would reach the
separatrix are independent of $C$ in this model. Using
$g_*=g_{*s}=61.75$ near the epoch gives $T_*=0.534$ GeV and
$H_*=3.05\times10^{-19}$ GeV. The prescribed classical seed also
uses angular variables independent of $f_a$. Once a physical
conversion exists, solve

$$
Y_\theta=\frac{d_a}{m_aC},\qquad
f_a^2=\frac{c_BT_{\rm EW}^2d_a}{m_aCY_B},\qquad
|g_{a\gamma\gamma}|=\frac{\alpha}{2\pi f_a}.
$$

These equations close the calibration algebra for that preparation.
A vacuum seed introduces the factor $(m_a/f_a)^2$ into its angular
power. A calculation including that seed must retain its dependence
when solving the matching. The linear vacuum-onset check varies the
normalization corresponding to $C=0.5,2,3$ as a sensitivity diagnostic.

The leading weak-anomaly logarithm gives a scale check. If the
fermion coefficients were set to zero instead at
$\Lambda=2.43\times10^8$ GeV, their left-handed quark and lepton
coefficients at the electroweak scale would be about
$-4.49\times10^{-4}$ by the cogenesis paper's leading-log formula.
This moves $c_B$ by about $-0.30\%$. The comparator retains the
zero boundary condition at 130 GeV. The on-shell weak-boson
and derivative-fermion loop terms have heavy-mass kinematic
suppression at 110 micro-eV, with $m_a^2/m_e^2=4.63\times10^{-20}$.
This is a size check, with full electroweak running, threshold
matching and the thermal sphaleron freeze-out profile left open.
[Co, Hall and Harigaya, Eqs. 2.8 and 3.2](https://arxiv.org/abs/2006.04809) ·
[Bauer et al., electroweak matching](https://arxiv.org/abs/2012.12272)

## Benchmarks and nonlinear evolution

The homogeneous proper-time ODE gives $C=2.00138$ at
$H_*/m_a=0.001$ and phase zero using the exact libration action.
Changing the initial phase to 1 and 2 radians in that slow
control keeps the result near 2. The rotating and librating
Lamé instability bands were checked against independent
monodromy integrations. At the rotating test point
$k/m_a=\sqrt{1.5}$, the rate is $0.207118m_a$ by monodromy and
$0.207124m_a$ from a weak wave in the 3D solver.
[Qualification data](results/qualification.json)

The radiation-era lattice evolves $\chi=a\theta$ in conformal
time $\tau=m_a\eta$. Here $a'=H_*/m_a$ and $a''=0$, so

$$
\chi''-\nabla^2\chi+a^3\sin(\chi/a)=0 .
$$

The projected sine force is evaluated on the periodic grid
with modes retained below the componentwise two-thirds cutoff.
The sine contains arbitrarily high harmonics, so momentum
headroom and spacing tests are essential. Matching Fourier
labels hold common random modes fixed between box resolutions.
The initial Gaussian state uses a specified fractional angular
velocity variance and an independent phase quadrature. Its
normalization is motivated by the size of primordial curvature,
while the transfer to the subhorizon ALP modes is a model input.

In units of $f_a^2m_a^2$, the simulation evaluates the full density
and pressure and checks $d(a^3\rho)/da=-3a^2P$. The late quadratic
dictionary uses the adiabatic action of $a^{3/2}\theta_k$ on a single
global vacuum branch. Its frequency squared in units of $m_a^2$
is $1+(k/a)^2+3H^2/(4m_a^2)$. The independent dictionary divides
the instantaneous quadratic free energy by the physical frequency.
Both omit nonlinear interactions. Their agreement leaves the
anharmonic energy and cutoff tail to independent checks. Number
estimation begins only in the registered late window.

| Run | Late quadratic proxy for $C$ | Proxy drift | Peak anharmonic energy | Peak cutoff tail | Energy residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| Base | 1.3294 | 2.66% | 5.57% | 11.12% | 0.887% |
| Spacing | 1.2401 | 1.99% | 5.26% | 4.70% | 0.986% |
| Volume | 1.3269 | 1.90% | 5.84% | 8.38% | 0.851% |
| Time | 1.3306 | 2.88% | 6.28% | 11.06% | 0.221% |
| Seed | 1.3333 | 2.64% | 6.03% | 10.93% | 0.885% |
| Lower seed | 1.4398 | 4.77% | 8.06% | 11.79% | 0.843% |
| Faster expansion | 2.0066 | 10.68% | 16.35% | 11.64% | 0.690% |
| Slower expansion | 1.0524 | 2.62% | 3.51% | 10.18% | 0.996% |

The base run transfers up to 76.4% of total energy into
gradient and inhomogeneous kinetic energy. The late zero-mode
fraction of the quadratic proxy is 0.40% at base and 0.56%
at finer spacing. The final momentum RMS is $0.428m_a$ and
$0.509m_a$ respectively in that dictionary. The complete
binned spectra and signed pressure work are archived. These
momentum statistics diagnose the simulation and do not set
a physical spectrum or halo line shape. The seed-amplitude
and hierarchy variations are sensitivities of unresolved proxies.
A physical error bar requires the physical hierarchy and preparation.
[Raw trace arrays and grid results](results/)

## Physical linear onset

The slowly swept Floquet calculation uses the exact rotating
Lamé rate along the adiabatic pendulum trajectory. It follows
the linear spectrum from $a=0.25$ to the first point where
predicted excess energy in fluctuations reaches one percent.
For the classical seed the peak physical momentum is
$7.16m_a$ there. For the canonical vacuum seeds it is about
$4.51m_a$. A direct linear integration across the principal
vacuum resonance gives log energy amplification 103.28,
compared with 101.93 in the swept approximation. This
cross-check establishes strong early growth and leaves the
subsequent nonlinear transfer open.

The dominant modes are well inside the horizon at electroweak
freeze-out, with $k/(aH)$ of order $3\times10^3$ to $4\times10^3$.
An adiabatic seed formula based on those modes being
superhorizon when charge generation starts needs a specified
earlier generation history. The fixed baryon and dark
densities alone do not determine that history.
[Fragmentation and initial condition domain](https://arxiv.org/html/2206.14259v2#S3.SS3) ·
[Lattice rescattering in the QCD axion comparator](https://arxiv.org/abs/2507.01822)

## Retrospective ORGAN comparison

After freezing the conversion decision, read the frequency-dependent
public tabulation of ORGAN Phase 1b Figure 5. The two mass rows
bracketing 110 micro-eV yield a nominal
$g_{95}\simeq4.5284\times10^{-12}\,\mathrm{GeV}^{-1}$
at 26.5979 GHz. Linear and logarithmic interpolation give the
same reported digits. This is a public curve reproduction,
with no bin-level likelihood reconstruction.

Under the paper's standard halo and all-local-dark-matter
assumptions, the leading EFT benchmark gives
$g(C)=4.7809\times10^{-12}\sqrt{C/2}\,\mathrm{GeV}^{-1}$.
The nominal crossing is **$C=1.794$**, compared with the
earlier 1.40 headline estimate. At $C=2$ the coupling is
5.58% above the reproduced curve. That margin is comparable
to the paper's reported 6.1% coupling systematic. This
is a retrospective comparison at a mass chosen within the
published window. It leaves the physical exclusion unresolved.
The local density and velocity distribution remain separate
inputs, and a prospective frequency interval is unopened.
[Experiment and systematic](https://arxiv.org/abs/2310.00904) ·
[Pinned curve data](https://github.com/cajohare/AxionLimits/blob/7d375f4879b32406a239fe48d2615a4bfd9bc0bb/limit_data/AxionPhoton/ORGAN-1b.txt)

A physical conversion calculation needs a specified preparation
through the linear instability, a matched nonlinear evolution
at the physical hierarchy, and a resolved late momentum tail
with a stable particle count. The first-percent growth times
set a concrete earlier matching target. The registered stop
condition applies before propagating any result to a new
experimental interval.

## Replay

From this directory use Python, NumPy and SciPy. The source
and gate scripts consume no detector data.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python qualify.py
python physical.py
for case in base spacing volume time seed lower_seed faster slower; do
    python lattice.py "$case"
done
python decide.py
# Optional independent curve replay after freezing the decision.
python organ.py /path/to/pinned/ORGAN-1b.txt
```

The NPZ ledger columns are $a$, comoving energy, its pressure
derivative, the gradient-plus-inhomogeneous-kinetic fraction,
and comoving mean angular velocity. The late arrays give
$a$, the two quadratic conversion proxies, zero-mode fraction,
momentum RMS, anharmonic fraction, branch failures, cutoff tail,
and $H/m_a$. These and the binned spectra are diagnostics where
the particle gate misses.

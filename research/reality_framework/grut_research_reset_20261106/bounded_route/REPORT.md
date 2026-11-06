# Physical ALP route closed and an independent vector benchmark retained

The bounded ALP test produced a converged linear covariance at
the physical hierarchy, then reached its cap without a stable
particle-counting regime. The next test excluded the frozen
40 GeV minimal Higgs-portal scalar as all of the dark matter from
equilibrium freeze-out. A separate, fixed 1 meV Proca vector
remains compatible with the BICEP/Keck tensor bound in the
prompt-radiation, nearly constant-inflation-Hubble limit.
Its required tensor-to-scalar ratio is about 0.0135.

[ALP protocol](PROTOCOL.md) · [ALP decision and raw results](results/route_decision.json) ·
[Scalar bound](results/portal.json) · [Vector protocol](VECTOR_PROTOCOL.md) ·
[Vector transfer](results/vector_transfer.json) · [Tensor comparison](results/vector_observation.json)

## The bounded ALP handoff

The frozen mass, anomaly texture and radiation model are the
ones in the [conversion packet](../alp_conversion/REPORT.md).
This test fixes a conditional canonical vacuum at $a=0.49$,
homogeneous phase zero, and $f_a=2.43\times10^8$ GeV for the
vacuum normalization. That value is an input from the
homogeneous conversion example. Physical conversion remains
undetermined. The physical $H_*/m_a$ is $2.774\times10^{-6}$.

Integrate the rotating mean and the two real fundamental
solutions of each retained perturbation,

$$
\chi_k''+(k^2+a^2\cos\bar\theta)\chi_k=0,\qquad
a'=H_*/m_a,\qquad \chi=a\delta\theta .
$$

Use the same two Hermitian Gaussian draws to assemble
$\chi_k$ and its canonical momentum. This preserves the
phase correlation of the amplified modes. The retained
band is $2\le k/m_a\le3$. The swept estimate puts earlier
amplified fluctuation energy at $3.24\times10^{-24}$ of
the homogeneous energy when the vacuum is prepared. The
same diagnostic puts omitted amplified energy below
$1.1\times10^{-29}$ of its estimate at both handoffs. The
exact finite-box covariance changes by at most
$4.87\times10^{-6}$ on tightening the ODE tolerance.
[Matched quadratures and volume diagnostics](results/match.json) ·
[Swept Lamé method](https://arxiv.org/abs/2206.14259)

At the two exact matched times $a=0.519471$ and $0.519770$,
the larger box samples the narrow shell differently.
Its expected fluctuation energy is 1.494 and 0.124 times the
reference box at the same times, with equal grid spacing.
Changing spacing at fixed volume preserves the initial
field and momentum on shared modes. The volume ratios
measure finite-box sampling of the linear state, with
its normalization held fixed.

| Nonlinear case | Endpoint status | Completed steps | Mean phase turns in the terminal window | Energy residual |
| --- | --- | ---: | ---: | ---: |
| Base | Reached | 30,000 | 66.25 | 0.0529% |
| Halved timestep | Reached | 60,000 | 66.44 | 0.0130% |
| Later handoff | Reached | 27,843 | 68.53 | 0.0444% |
| Finer spacing | Wall cap | 28,590 of 30,000 | Partial | 0.0519% |
| Larger volume | Wall cap | 28,444 of 30,000 | Partial | 0.0504% |

The three complete runs retain 75.6–81.4% of their final
energy in the mean field's kinetic term. The base reaches
a gradient and inhomogeneous kinetic share of 21.9%.
The mean still winds through the final 100 conformal
units, and the quadratic proxy drifts by 11–13% in the
complete runs. All three fail the fixed-branch test. The
larger runs preserve partial histories through the cap.
Only common endpoints enter the terminal proxy comparison.
Energy closure and a successful linear match therefore
leave the late particle conversion open. The resource
rule closes this lane at the registered window and leaves
the prospective comparison unopened.
[Prescribed gate and all measured misses](results/route_decision.json)

## A halo-independent exclusion at 40 GeV

Use a stable real scalar with the single Higgs portal
$V\supset\lambda S^2H^\dagger H/2$ and conventional
annihilation into SM states during equilibrium freeze-out.
At $m_S=40$ GeV the portal gives

$$
\Gamma(h\to SS)=
\frac{\lambda^2 v^2}{32\pi m_h}\sqrt{1-\frac{4m_S^2}{m_h^2}} .
$$

The published ATLAS invisible branching limit of 0.107
at 95% confidence gives $|\lambda|\le0.01151$ when the
SM width is 4.10 MeV and the full invisible allowance
is assigned to $SS$.
[ATLAS invisible Higgs result](https://arxiv.org/abs/2301.10731) ·
[SM width input](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/9da32dcda3b65d4763b4af191d553018247adcfa/data/lhc-hxswg/sm/br/BR.txt)

A bound on the remaining depletion suffices. At $x=m_S/T=20$,
the abundance from an equilibrium start is at least the
equilibrium yield computed with the upper entropy count
106.75. For all earlier $x\le20$ that lower bound lies
below the actual equilibrium yield. For later temperatures
the Boltzmann equation implies

$$
\frac{1}{Y_\infty}\le
\frac{1}{Y_{\rm eq,min}(20)}
+\int_0^{2\,\mathrm{GeV}}
  \frac{s\,\langle\sigma v\rangle}{HT}
  \left(1+\frac{T}{3}\frac{d\ln g_{*s}}{dT}\right)dT .
$$

The integrand takes the maximal annihilation rate after dropping
inverse production. Use radiation domination, conserved
entropy, $g_{*s}/\sqrt{g_{*\rho}}\le12$ and an entropy-change
integral bounded by $T_0\ln(106.75/3.90)/3$. The annihilation
kernel is the inclusive virtual-Higgs rate

$$
\sigma v_{\rm cm}(E)=
\frac{2\lambda^2v^2\Gamma_{\rm SM}(E)}
 {E[(E^2-m_h^2)^2+m_h^2\Gamma_h^2]} .
$$

Inflate the virtual Higgs width to 20 MeV from 80 to 140
GeV, where the pinned table peaks at 8.17 MeV. Near the
pole use a lower propagator width of 4.0 MeV and discard
it away from the pole. The relativistic Maxwell–Boltzmann
average uses the exact single-integral kernel. Below 95 GeV
use the largest allowed cross section in that range.
Above 95 GeV the Bessel bounds
$K_2(x)\ge\sqrt{\pi/(2x)}e^{-x}$ and
$K_1(y)\le\sqrt{\pi/(2y)}e^{-y}(1+3/(8y))$ give a
weight upper bound that decreases for $x\ge20$.
Keep the pole in the integration. Above 140 GeV use a
loose $\sigma v_{\rm cm}$ ceiling of $1\,\mathrm{GeV}^{-2}$.
The s-channel partial-wave bound is far below that
ceiling there, as is perturbative $SS\to hh$ at the
tested couplings. This last tail adds
$2.93\times10^{-12}\,\mathrm{GeV}^{-2}$ to the
rate bound.
[Scalar singlet annihilation framework](https://doi.org/10.1140/epjc/s10052-017-5113-1) ·
[Relativistic thermal average](https://doi.org/10.1016/0550-3213(91)90438-4)

At the collider ceiling the uniform bound is
$\langle\sigma v\rangle\le9.49\times10^{-11}\,\mathrm{GeV}^{-2}$.
Using the tabulated widths gives $4.95\times10^{-12}$
at $x=20$, a check on the margin. The late yield is at
least $5.23\times10^{-11}$, or **4.76 times** the target
$1.10\times10^{-11}$. Inverting the conservative bound
requires an invisible branching fraction of at least
**0.410** for all the observed dark matter. This conflicts
with ATLAS at the fixed mass within the standard thermal
history and minimal portal. There is no halo-density input.
[Cosmological abundance](https://arxiv.org/abs/1807.06209)

## A surviving tensor test for a fixed vector

The separate analytic screen fixes a stable decoupled
spectator Proca mass of 1 meV. A relativistic longitudinal
mode is generated by the inflationary vacuum, with
initial $A_L$ power $(kH_k/(2\pi m))^2$ in the
slowly varying inflation limit. Use prompt radiation after
inflation and standard entropy conservation. No adjustable
interaction sets the yield.
[Inflationary production mechanism](https://arxiv.org/abs/1504.02102)

With $q=k/(a_*m)$ and $H(a_*)=m$, the radiation equation for
$u=A_L/A_{L,\rm in}$ in $x=a/a_*$ is

$$
u_{xx}+\frac{2q^2}{x(q^2+x^2)}u_x+(q^2+x^2)u=0 .
$$

Set $b=x/\sqrt{q^2+x^2}$, $z=bu$, and
$\omega^2=q^2+x^2-b_{xx}/b$. The late canonical
action gives $W(q)=q^2(z_x^2+\omega^2z^2)/\omega$
per logarithmic momentum. The infrared coefficient
$W/q^2=2\Gamma(3/4)^2/\pi$ follows from the zero-momentum
Bessel solution. In the relativistic ultraviolet,
$u\simeq\sin(qx)/(qx)$ gives $W\simeq1/q$.
The fine calculation agrees with the coefficients to
$2.58\times10^{-4}$ and $9.5\times10^{-7}$ relative.
The integrated weight is $I=1.18142$. It changes by
$6.73\times10^{-5}$ with the grid, cutoff and initial-time
refinement, and $1.23\times10^{-4}$ between late readouts.
The corresponding scale is $T_*=844$ GeV and
$k_*=1.45\times10^{10}\,\mathrm{Mpc}^{-1}$.

The energy per radiation entropy is
$m^2 H_{\rm eff}^2 I/(8\pi^2 s_*)$ with $H_{\rm eff}^2$
the average of $H_k^2$ weighted by $W(q)$. The observed
dark matter density requires $H_{\rm eff}=2.88\times10^{13}$
GeV. For nearly constant inflationary $H$, ordinary
vacuum tensors give $r_{0.05}=0.01345$ at scalar power
$2.1\times10^{-9}$. The production result was frozen
before comparing with the published BK18 limit
$r_{0.05}<0.036$ at 95% confidence.
[Published production normalization, Eq. 46](https://arxiv.org/pdf/1504.02102) ·
[Vacuum tensor normalization](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-inflation.pdf) ·
[BICEP/Keck XIII](https://arxiv.org/abs/2110.00483)

The modes that dominate production exit much later
than the CMB pivot, separated by about 26 in log
wavenumber. For a slowly varying, nonincreasing $H$,
the pivot ratio is at least approximately 0.01345
if this fixed mass makes all the dark matter. Earlier
modes are negligible in the transfer for sub-Planckian
inflation. BK18 requires $H_{\rm eff}/H_{\rm pivot}>0.611$
at this benchmark. Its corresponding lower mass is
approximately $0.140$ meV for this radiation history
and constant relativistic degrees of freedom. The
mass is a comparator choice, and the CMB comparison
provides a test of its compatibility. A primordial
tensor upper limit below the required amplitude would
exclude this fixed point under those production assumptions.

## Replay

Use Python, NumPy and SciPy. The wall caps in the ALP
protocol are part of the test. The archived partial
histories retain their actual stopping steps.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python match.py
printf '%s\n' base spacing volume time later | \
    xargs -P 5 -I CASE timeout --kill-after=15s 360s python evolve.py CASE
python decide.py
python portal.py /path/to/pinned/BR.txt
python vector_screen.py
```

The matched NPZ contains the mode labels, times and
background-plus-quadrature states. Nonlinear ledgers
store conformal time from the first seam, scale factor,
comoving energy, its pressure derivative, gradient and
nonzero-mode kinetic fraction, coherent kinetic fraction,
and mean phase. The vector spectrum is $W(q)$ on the fine
momentum grid. The observational file uses the frozen
$r_{\rm flat}$ and BK18 upper bound, with
$m_{\min}=m(r_{\rm flat}/0.036)^2$ and
$H_{\rm eff}/H_{\rm pivot}>\sqrt{r_{\rm flat}/0.036}$.

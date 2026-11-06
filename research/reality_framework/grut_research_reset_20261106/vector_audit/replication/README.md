# Fixed 1 meV vector replication packet

The conditional counterexample survives a separate implementation at the
frozen parameters. At the fitted spectra the smooth stiff-to-prompt
abundance ratios agree with the archive within 0.01%. With a **common
primordial tensor spectrum at the same comoving wavenumber**, the stiff
history gives **4.950 times** the early-radiation abundance. This matched
comparison changes only the transfer applied to that input spectrum.
The case remains conditional on the specified expansion and conversion.

| Smooth stiff case | Archived fitted ratio to prompt | Separate calculation | Relative difference | Ratio with common stiff spectrum |
|---|---:|---:|---:|---:|
| Narrow blend | 0.999972 | 0.999880 | \(-0.00917\%\) | 4.94985 |
| Wide blend | 0.999966 | 0.999880 | \(-0.00855\%\) | 4.94985 |

The matter-era control followed by early radiation gives the same
matched comparison to the shown digits. The flat-spectrum comparison
gives 4.961. The fixed fitted ratios near unity reflect the separate
original abundance calibrations. The common-spectrum comparison
isolates the reheating dependence. The original stiff bracket
\(\alpha\in[1,40]\) missed the abundance target. The frozen
\([0.5,1]\) extension gave \(r_{0.05}=0.002765\), with the two
early-radiation fits near 0.0148. No parameter was refitted here.
The instantaneous-transition failure is retained in the original audit.
[Original audit](https://github.com/ryangrvr/TestingGRUT/blob/97e1772c126d05309d2f8dc36ec65780de275b19/research/reality_framework/grut_research_reset_20261106/vector_audit/CONSISTENCY.md).

## Implementation and normalization

The inputs preserve the canonical E-model potential, minimal Abelian
Proca field, growing vacuum mode, measured scalar calibration used by
the original fit, ideal thermal degrees of freedom and two smooth
background prescriptions. The code uses the same physical equations
and general NumPy and SciPy tools. It imports no archive solver.
It integrates the scalar and Friedmann equations in cosmic time,
integrates the smooth post-inflation background from conservation,
and evolves canonical inflationary modes and the second-order
longitudinal equation in the linear scale factor. The archive used
log scale factor and a first-order vector Hamiltonian system.
The new readout averages physical mode energy divided by frequency
in two late windows.

Eliminating the constrained Proca time component gives one longitudinal
mode with kinetic coefficient \(F=a^2m^2/(k^2+a^2m^2)\). At relativistic
inflationary exit \(m A_L/k\) is the massless spectator amplitude
\(\varphi\). Two tensor helicities give \(P_t=8P_\varphi/\bar M^2\)
with \(\bar M\) the reduced Planck mass. This relation uses the same
canonical vacuum normalization per degree of freedom. The two transverse
Proca states are conformal in the relativistic limit. Their leading
inflationary amplification vanishes in that limit. [Production
mechanism](https://arxiv.org/abs/1504.02102).

For \(x=a/a_\ast,\ q=k/k_\ast,\ h=H/m\), let \(L\) denote the
longitudinal transfer of a unit growing mode and \(D=x^2h\,dL/dx\).
The physical energy and isotropic pressure have common factor
\(1/(2a^4)\) in conformal units, with quadratic numerators
\(E=F D^2+x^2L^2\) and \(p=[(3-2F)FD^2-x^2L^2]/3\).
The mode equation gives \(xE_x-E+3p=0\), recovered numerically
to \(4\times10^{-12}\) relative on the resolved vector grid.
With \(\omega=\sqrt{x^2+q^2}\), the late averaged \(W=q^2E/\omega\)
gives \(d_A=m^2\bar M^2\int W P_t\,d\ln q/(16 S_\ast)\).
Here \(S_\ast=s_0/a_\ast^3\) in the convention \(a_0=1\).
The independent pivot solution \(k=aH\) gives 60.6424 and 60.6252
e-folds to the smoothed end. Entropy and \(k_\ast=a_\ast m\) agree
with the archive. The tensor normalization at the frozen scalar
amplitude agrees to \(2.1\times10^{-8}\) relative. These agreements
measure numerical consistency with shared inputs.

## Spectrum and stress audit

The two grids evolve the soft longitudinal modes directly. The Bessel
continuation from \(q=128\) to the archived upper cut contributes 11.7%
of the stiff abundance. Its coefficient agrees with the direct solution
at that edge to \(3\times10^{-6}\) relative. Beyond \(q=10^8\), the
growing-mode tail is bounded by \(3.4\times10^{-4}\) of the target
under the nonincreasing-power continuation. The independently
bounded hard transition modes can contribute under \(9.1\times10^{-6}\)
more in rest energy. The fitted ratio mesh shifts and late-window shifts
are recorded separately from these tail estimates.

For the hard modes we match the independently evolved canonical solution
to the real stiff Bessel basis. On \(0.001\le K=k/k_{\rm end}\le64\),
the resolved \(K^4|\beta_K|^2\) energy integrals agree with the archive
to \(9\times10^{-7}\) relative. For the omitted band use the instantaneous
frequency \(\omega_K=\sqrt{K^2-U}\), \(U=a''/a\) in units of
\(k_{\rm end}^2\). The exact mixing coefficient is
\(\omega_K'/(2\omega_K)\). Integrate its leading phase integral
twice by parts and bound the nonlinear remainder by
\(J(\cosh J-1)\), with \(J=\int|\omega_K'/2\omega_K|d\eta\).
The derivative quadrature, tested at two resolutions with a
factor-two margin, bounds the additional hard energy by 2.86% and
0.0831% of the resolved integrals. Both are below the registered 5%.
For a sharp step the leading \(K^4|\beta_K|^2\) instead approaches
a positive constant, which gives the logarithmic divergence.

One relativistic longitudinal degree of freedom and two tensor
helicities yield a combined conservative \(\Delta N_{\rm eff}\)
upper estimate of 0.00111 for the two smooth cases. This estimate
uses the shared low-momentum gap envelope and the same thermal entropy
normalization. The new stress check covers the longitudinal stress
identity and the produced hard energy after oscillation averaging.
The archive's time-resolved hard vacuum subtraction and transient
pressure envelopes remain separate. BBN Table 2 with PDG abundances
gives \(-0.10\pm0.21\) for PArthENoPE v3.0 and
\(-0.14\pm0.21\) for PRyMordial with NACRE II.
[BBN Table 2](https://arxiv.org/html/2401.15054v3).

## Perturbative scope and main physical failure mode

One effective closure uses two barotropic perfect fluids with
\(w_s=3/5,\ w_r=1/3\) and energy-momentum exchange \(\pm Q u^\mu_{\rm tot}\).
Choose the local scalar
\(Q=(\sigma-0.8)\sqrt{\rho/(3\bar M^2)}\,\rho_s\rho_r/\rho\), where
\(\rho=\rho_s+\rho_r\). It reproduces the frozen logistic conversion
in a homogeneous solution and conserves total energy-momentum.
At linear order
\(\delta Q/Q=\delta\rho_s/\rho_s+\delta\rho_r/\rho_r-\delta\rho/(2\rho)\).
A common time shift gives \(\delta Q=\dot Q\,\delta t\) and
\(\delta p_{\rm nad}=0\). A separate-universe integration checks
the adiabatic mode on uniform-density surfaces before, during
and after conversion. This illustrates a local closure that
admits the assumed long-wavelength adiabatic solution.
[Adiabatic energy transfer](https://arxiv.org/abs/astro-ph/0211602).

The main physical vulnerability is an independent transfer or radiation
seed perturbation. For this example, with \(f=\rho_r/\rho\), a seed
change at fixed initial density gives
\(\delta N\simeq (0.2/\sigma)\,\delta\log[f/(1-f)]\) on a late
uniform-density surface. The numerical test recovers the response.
This closure requires a fractional radiation seed of order
\(10^{-260}\) or \(10^{-130}\) at the end blend for the two rates.
These tiny seeds illustrate the dependence of this example on the
initial bath. The homogeneous transfer leaves the intrinsic transfer
perturbation and the physical origin of the bath open. A different
physical completion can alter the scalar perturbation and therefore
the interpretation of the fitted \(r\). The replication is by a
separate implementation in this packet. External reproduction
and any CMB tensor or tilt holdout remain separate tasks.

## Reproduce without the original repository

Use Python 3.12 and the tested dependencies below. The measured arrays,
comparison values and protocol are included. Run in the extracted
packet directory.

```sh
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
for history in prompt early_matter narrow wide; do
  python replicate.py "$history"
  python replicate.py "$history" --fine
done
python check.py
```

The runs write [measured outputs](measured/). The final command checks
[the frozen references](archived_reference.json) after the numerical
measurements are present. Read [the protocol](PROTOCOL.md),
[the comparison](comparison.json), and the [separate audit checks](independent_checks.json).

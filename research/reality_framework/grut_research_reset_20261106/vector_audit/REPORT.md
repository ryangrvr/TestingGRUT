# When a 1 meV inflationary vector requires observable CMB tensors

## Verdict

The all-dark-matter requirement connects the vector abundance to a positive
weighted integral of the primordial tensor spectrum. For the stable 1 meV
minimal Proca vacuum benchmark, radiation transfer and a nonincreasing tensor
spectrum after the CMB pivot give \(r_{0.05}\gtrsim0.01345\). The two
calibrated histories with radiation in place before the relevant modes evolve
give \(r_{0.05}=0.01473\) and \(0.01486\).

A post-inflationary fluid with \(w=3/5\) that reheats at 200 GeV invalidates
the radiation floor. Even the resolved \(0.25\le k/k_\ast\le 4\) band produces
the full required density at \(r=0.006074\) in the specified E-model family.
The late-stiff calculation tests a distinct reheating history at the same
mass. The tensor relation is a conditional prediction of the published
gravitational-production mechanism. No entry in the GRUT theory tuple
selects the vector mass, inflation potential, or reheating history.

## Normalization and abundance

Use reduced Planck mass \(M=2.435\times10^{18}\) GeV and scale factor \(a_0=1\).
Start with \(-F_{\mu\nu}F^{\mu\nu}/4-m^2 A_\mu A^\mu/2\).
In conformal time \(\eta\), the auxiliary field is
\(A_0=-ik A_L'/(k^2+a^2m^2)\) with Fourier convention \(e^{i\mathbf{k}\cdot
\mathbf{x}}\) and \(A_L=\hat{\mathbf{k}}\cdot\mathbf{A}\). Eliminating it gives

\[
 S_L=\frac12\int d\eta\,\frac{d^3k}{(2\pi)^3}
 K\left(|A_L'|^2-(k^2+a^2m^2)|A_L|^2\right),
 \qquad K=\frac{a^2m^2}{k^2+a^2m^2}.
\]

The canonical field is \(B=\sqrt K A_L\). For \(k/a\gg m\) the field
\(\varphi=B/a\simeq(m/k)A_L\) obeys the massless spectator equation.
The vacuum spectrum in that regime is
\(\mathcal P_{A_L}=(k/m)^2\mathcal P_\varphi\).
This counts one longitudinal polarization and integrates each Fourier mode
once. The other two Proca polarizations approach the conformal Maxwell
vacuum for \(m/H\to0\). They supply no leading inflationary growing spectrum.
In Einstein gravity an orthonormal helicity tensor gives the canonical
variable \(aM h_\lambda/2\) for each metric amplitude \(h_\lambda\).
Summing the two helicities gives
\(\mathcal P_t=8\mathcal P_\varphi/M^2\). This proportionality uses the
mode equation for any of the backgrounds tested here, without setting the
Hubble rate equal at different exit times.

Set \(H(a_\ast)=m\), \(x=a/a_\ast\), \(h=H/m\), \(q=k/(a_\ast m)\) and
\(n=\ln x\). For unit growing mode \(X\) and its dimensionless momentum
\(P=xKh\,dX/dn\), the transfer equations are
\[
 \frac{dX}{dn}=\frac{P}{xKh},\qquad
 \frac{dP}{dn}=-\frac{xX}{h}.
\]
Continuity of \(X\) and \(P\) matches the solutions through reheating.
Compute the late quadratic action of \(B_X=\sqrt{K}X\) per mode in
\(\hat\eta=a_\ast m\eta\). Write the dimensionless integrand
\(W(q)=q^2(|dB_X/d\hat\eta|^2+\hat\omega^2|B_X|^2)/\hat\omega\), with
\(\hat\omega^2=q^2+x^2-(\sqrt K)_{,\hat\eta\hat\eta}/\sqrt K\).
This counts twice the action before multiplication by \(q^2\). We evaluate it
in the late adiabatic regime and use the massive rest energy for the
late-time density. The frozen classical modes have high occupation and
the unamplified vacuum contribution is omitted.

Let \(S_\ast=s_R(a_R/a_\ast)^3\) be the radiation entropy per \(a_\ast^3\).
It equals \(s(T_\ast)\) when \(H=m\) occurs in adiabatic radiation. Then
\(a_\ast=(s_0/S_\ast)^{1/3}\) and \(k_\ast=a_\ast m\) in natural units.
The abundance and tensor relation are
\[
 d_A\equiv\frac{\rho_A}{s}
   =\frac{m^2 M^2}{16S_\ast}\int_0^\infty
       W(q)\mathcal P_t(k_\ast q)\,d\ln q,\qquad
 r_{0.05} F = r_{\rm flat}
   =\frac{16d_A S_\ast}{m^2M^2 A_s I},
\]
\[
 I=\int W\,d\ln q,\qquad
 F=\frac{1}{I}\int W(q)
      \frac{\mathcal P_t(k_\ast q)}{\mathcal P_t(k_p)}\,d\ln q .
\]
These identities apply at leading order in the relativistic-vacuum,
linear spectator and late adiabatic approximations. The relation
\(F\le1\) follows if \(\mathcal P_t(k)\le\mathcal P_t(k_p)\) throughout
the contributing range after the pivot. Monotonic \(H\) supports that
condition in smooth vacuum slow roll. The bound here uses the tensor
spectrum condition explicitly. Earlier modes have \(4.84\times10^{-24}\)
of the radiation transfer weight, for a possible absolute \(r\) correction
of only \(2.3\times10^{-15}\) if their power stays below unity.

The radiation equations give \(T_\ast=843.61\) GeV, \(S_\ast=2.81135\times
10^{10}\ {\rm GeV}^3\), and \(k_\ast=1.4458\times10^{10}\ {\rm Mpc}^{-1}\)
at \(g_\rho=g_s=106.75\). Using \(T_{\gamma0}=2.7255\) K, \(g_{s0}=3.91\)
and \(\rho_c/h_0^2=3M^2H_{100}^2\), the conversion
\((\Omega_c h_0^2)\rho_c/(h_0^2 s_0)\) is 0.4371 eV at
\(\Omega_c h_0^2=0.120\). The frozen 0.44 eV input rounds this by 0.66%.
In radiation \(S_\ast\) scales as \(g_s g_\rho^{-3/4}\). The ideal
thermal bath, later entropy conservation, and all-vector dark matter are
assumptions of the quoted floor.

## Independent checks of transfer and inflation

In radiation at \(q\ll1\), \(X\propto\sqrt x J_{-1/4}(x^2/2)\)
with regular initial data. This gives
\(\lim_{q\to0} W/q^2=2\Gamma(3/4)^2/\pi\).
At high \(q\), in a phase with \(h=x^{-p}\), \(p=3(1+w)/2\),
the relativistic solution is
\(2^\nu\Gamma(1+\nu)y^{-\nu}J_\nu(y)\) for
\(y=qx^{p-1}/(p-1)\) and \(\nu=(3-p)/(2(p-1))\).
Its adiabatic action gives
\(W\simeq [2^{2\nu+1}\Gamma(1+\nu)^2(p-1)^{2\nu+1}/\pi]q^{-2\nu}\).
Here \(2\nu=1\) for radiation and \(3/7\) for the late-stiff case.
The low-\(q\) radiation coefficient agrees to \(2.58\times10^{-4}\).
The high-\(q\) coefficient agrees to better than \(1.4\times10^{-5}\)
at \(q=128\) in the three controls.

The radiation integral is 1.18142414. An independently implemented
conformal-time transfer gave 1.18142416. This agreement is finer than
the \(6.74\times10^{-5}\) grid shift and \(1.23\times10^{-4}\) late readout
shift, which are the relevant numerical checks. The stiff integral is
2.85496 with grid and readout shifts \(2.48\times10^{-5}\) and
\(1.72\times10^{-5}\). Its \(q>128\) asymptotic continuation carries
11.8% of the flat-power integral. The finite-band counterexample removes
that continuation entirely. The first evolving-inflation attempt at
\(q_{\rm stop}=10^8\) hit the registered minimum freeze-out interval for
a matter-era bracket endpoint. The replay ends the tensor integration
at \(q=10^6\) for the radiation tail and \(10^8\) for the stiff tail,
or earlier to leave at least eight e-folds after exit, and records the
residual weight separately. This changes no fitted input or gate.

For the E-model runs we solve the canonical background to \(\epsilon_H=1\),
match \(a_{\rm end}\) using the entropy and the declared reheating history,
and derive the pivot exit from \(k_p=aH\). The scalar Mukhanov–Sasaki
and tensor vacuum modes are integrated on this background. Independent
constant-\(\epsilon\) Hankel solutions give relative errors below
\(3\times10^{-9}\), including the power-law scalar control. Moving the mode
start, the late readout and the quadrature nodes tests numerical sensitivity.
The amplitude \(A_s=2.1\times10^{-9}\) and \(d_A=0.44\) eV set the two
E-model parameters for each abundance fit. The tilt, its running and the
tensor ratio are outputs.

| History | \(\alpha\) | Pivot exit e-folds before end | Representative vector exit e-folds before end | \(F\) | \(r_{0.05}\) | \(n_s\) |
|---|---:|---:|---:|---:|---:|---:|
| Prompt radiation | 5.78981 | 56.31 | 29.56 | 0.91361 | 0.0147266 | 0.96727 |
| \(w=0\) to \(10^9\) GeV | 4.73259 | 51.36 | 24.60 | 0.90517 | 0.0148640 | 0.96403 |
| \(w=3/5\) to 200 GeV | 0.93340 | 60.66 | 33.70 | 0.98051 | 0.0027648 | 0.96804 |

The representative vector mode uses \(q=1.375\). The spectrum weights a
range around it. In the radiation cases \(k_\ast/k_p=2.89\times10^{11}\),
so the representative mode exits about 27 e-folds after the pivot.
Using the pivot tensor power at all of these scales would underestimate
the required \(r\) by 9.5% and 10.5% relative to the flat-power value.
The \(A_s\)-only controls with prompt radiation produce fractions 0.2515
at \(\alpha=1\) and 0.7721 at \(\alpha=4\), at \(r=0.003463\) and
0.01114. These controls check the response to changing the inflation
history without calibrating the abundance. The two radiation fits
illustrate specified histories within the broader lower-bound class.

The original late-stiff bracket \([1,40]\) missed because the endpoint
densities were 1.063 and 11.97 times the target. Its finite-band fit
bracketed with \(r=0.006074\) and a full calculated density 2.15 times
the target at that point. The separately registered \([0.5,1]\)
extension kept the same mass, background and matching, and bracketed
with endpoint densities 0.568 and 1.063. Its fine fit gives the third
row of the table. The direct finite band contributes 46.4% there and the
unincluded far tail is below \(3.4\times10^{-4}\) of the target under the
decreasing-power continuation. This tail estimate limits the precision
of the abundance calibration. Numerical refinement of all three
bracketed fits shifts \(r\) by less than \(4.2\times10^{-8}\) relative
and \(n_s\) by less than \(2\times10^{-10}\) absolute. The transfer
checks and modeling approximations set the larger uncertainty.

## When the reheating assumption matters

The transfer requires a sufficiently long relativistic,
superhorizon, frozen interval before horizon entry. Here \(q_R=a_R
H_R/(a_\ast m)\) measures where modes meet the reheating transition.
Prompt radiation starts at \(q_R\simeq3.0\times10^{12}\). In the early
matter case \(H_R/m=1.405\times10^{12}\) and \(q_R=1.185\times10^6\).
For the modes with \(q\ll q_R\), gradient corrections at that transition
scale as \((q/q_R)^2\). The radiation tail fraction above \(q_R\)
is about \(7.1\times10^{-7}\) in the early matter case. Those earlier
entrants have a steeper high-\(q\) transfer under \(w=0\).
The two radiation transfer integrals agree to \(2.73\times10^{-11}\)
on the computed range with the same radiation-tail continuation.
The high-\(q\) weight estimate sets the scale of the omitted earlier
entry correction. Early matter changes the pivot matching even
though the dominant transfer remains radiation-equivalent.

For \(w=3/5\) to 200 GeV, \(H_R/m=0.0562\), \(q_R=0.1865\), and
\(S_\ast=1.36886\times10^{10}\ {\rm GeV}^3\). The main production
modes evolve before the radiation bath exists. The \(W\) integral is
2.42 times its radiation value and \(S_\ast\) is 0.487 times its
radiation value. Both effects change the density per unit tensor power.
By the transition to radiation the smallest oscillation frequency is
about \(18 H_R\). The flat-power tensor requirement is 0.002711, and
the evolving-inflation result is 0.002765. The finite-band result shows
that the violation of the radiation floor survives removal of the
extrapolated ultraviolet modes.

The stiff period spans about 25 e-folds. As an explicit effective
history it uses a conserved dominant component with pressure
\(3\rho/5\), followed by energy conversion on a uniform-density surface
with adiabatic perturbations. A separate canonical scalar on the
scalar-dominated exponential scaling solution with slope
\(\sqrt{24/5}/M\) can supply that equation of state. The E-model
inflaton has a quadratic minimum, so extending that potential alone
through reheating would give the matter-like case. Rapid transitions,
the conversion time, adiabatic matching and absence of additional
hard production are supplied assumptions. The calculation leaves the
microphysics of the conversion and detailed thermalization open.

## Spectator and stability domain

The minimal vector has positive \(m^2\) and \(K>0\) for \(a>0\).
A free Abelian Stückelberg description can preserve its constant mass
and its exact vector sign symmetry. Assume an EFT cutoff above the
inflationary scales and omit direct inflaton, curvature and
visible-sector portals. With the vector as the only odd state, the
symmetry forbids a single-vector decay. At the fitted scales
\(m/H\lesssim8\times10^{-26}\). The plateau estimate
\(H^4/(3M^2H^2)\) stays below \(5.8\times10^{-11}\) during inflation.
This compares a stress of order \(H^4\) with the background, with the
unamplified ultraviolet vacuum renormalized. For example, the leading
nearly de Sitter superhorizon mass contribution integrates to order
\(H^4/(16\pi^2)\). The modes that dominate the abundance also have
\(k/a\gg m\) at inflation end.

The initial stiff calculation used a horizon-cutoff stress estimate.
The consistency audit below supersedes that estimate and resolves
the transition-generated tensor and vector energy. The laboratory
interaction strength is unspecified beyond the benchmark's
minimal gravitational coupling.

## Stiff-history consistency

The original sharp history has a cutoff-dependent produced stress.
Two specified smooth, energy-conserving backgrounds preserve the
frozen abundance at the same \(r\), with small amplified vector and
gravitational-wave stress. The separate consistency audit provides
the UV, extra-radiation, pivot and cold-matter checks. The physical
counterexample remains conditional on the effective reheating
prescription. The previously seen CMB tensor and tilt marginals
cannot be held out. [Consistency audit](CONSISTENCY.md).

The Planck base-model uncertainties of 1% in \(\Omega_c h_0^2\)
and 1.4% in \(A_s\) enter the floor as \(d_A/A_s\). Ignoring covariance
gives a 1.7% calibration scale, separate from thermal, inflation and
reheating choices. Rescaling only to the Planck central values gives
a radiation floor of 0.01337. The inflation correction in the two
radiation fits is about 10%. The late-stiff history moves \(r\) by a
factor of about 5.3 relative to the prompt E-model fit.
[Planck Table 5](https://www.aanda.org/articles/aa/pdf/2020/09/aa33910-18.pdf).

A future tensor upper limit sufficiently below \(0.0134\) would
exclude the all-dark-matter benchmark within the stated radiation
transfer and decreasing-power class. The stiff counterexample
shows that this test depends on reheating. A tensor detection
at the calculated values would remain compatible with many mechanisms.

## Reproduction and program status

From this directory with Python, NumPy and SciPy, run

```sh
python transfer.py
OPENBLAS_NUM_THREADS=1 python inflation.py
OPENBLAS_NUM_THREADS=1 python inflation.py --late-extension
python diagnostics.py
```

The transfer and all three saved spectra replayed byte for byte in
a separate directory. The original registered bracket miss and the
extension are archived separately. [Protocol](PROTOCOL.md).
[Transfer](results/transfer.json). [Inflation](results/inflation.json).
[Extension](results/late_stiff_extension.json).
[Scale diagnostics](results/diagnostics.json).

Retain the 1 meV vector as a conditional comparator within standard
gravity and quantum field theory. The ALP conversion route remains
closed as computationally unresolved. The 40 GeV minimal thermal
scalar result remains scoped to the stated equilibrium freeze-out
history and all-dark-matter assumption, pending independent review.
[Earlier benchmark report](../bounded_route/REPORT.md).

The production baseline is Graham, Mardon and Rajendran
[arXiv:1504.02102](https://arxiv.org/abs/1504.02102), with reheating
dependence of gravitational production studied by Kolb and Long
[arXiv:2009.03828](https://arxiv.org/abs/2009.03828). The inflation
test uses the E-model \(\alpha\)-attractor potential specified in
the protocol. [Carrasco et al.](https://arxiv.org/abs/1504.05557).

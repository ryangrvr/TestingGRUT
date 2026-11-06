# Physical consistency of the 1 meV stiff-history counterexample

## Decision

The frozen 1 meV stiff fit has \(r_{0.05}=0.00276479\). The instantaneous
history fails the renormalized stress test. Its two jumps in the conformal
mode frequency produce a logarithmically divergent excitation energy.
Two energy-conserving smooth effective backgrounds preserve the frozen
abundance and clear the registered spectator, extra-radiation and cold
matter checks. This establishes a conditional counterexample to the
early-radiation tensor floor. A microscopic model for the conversion
to visible radiation and a distinct observational test remain open.

The mass, inflation potential parameters, scalar amplitude calibration
and radiation entropy are the same throughout. The original \([1,40]\)
abundance bracket miss and the separate \([0.5,1]\) extension remain
archived. The extension used \(\alpha=0.9333995236\) and
\(V_0/M^4=8.789272149\times10^{-11}\).
[Registered consistency test](CONSISTENCY_PROTOCOL.md).
[Frozen fit and bracket](results/late_stiff_extension.json).

## The sharp background fails at high momentum

For each tensor helicity and the relativistic longitudinal vector,
the canonical mode satisfies \(v''+(k^2-U)v=0\) with \(U=a''/a\).
Continuity of \(a,H,v,v'\) across a pressure step permits a jump in
\(U=(aH)^2(2-\epsilon_H)\). At high \(k\) the resulting
\(|\beta_k|=|\Delta U|/(4k^2)+O(k^{-3})\). The produced-energy
integral contains \(k^4|\beta_k|^2\,d\ln k\). Two distinct steps
contribute a positive sum of squares to the logarithmic divergence.
The cross term oscillates. The divergence is in the excited state,
so subtracting the unamplified vacuum leaves it in place.
A cutoff at the last inflationary horizon exit would conceal it.

## Matched smooth histories

Use \(N=\ln(a/a_{\rm end}^{\rm sharp})\) and keep the original
inflation solution before the first blend. Across \(N=\pm\delta\)
blend \(\ln H\) between the continued canonical inflation solution
and the stiff solution with a regularized beta polynomial of
degree nine. The first four blend derivatives vanish at each edge.
The smoothed \(\ln H\) joins both original trajectories exactly.
During this blend \(0<\epsilon_H<3\). A canonical scalar along the
trajectory realizes it with \(|d\phi/dN|=M\sqrt{2\epsilon_H}\) and
\(V=M^2H^2(3-\epsilon_H)\). The stiff segment corresponds to the
exponential scaling solution with slope \(\sqrt{24/5}/M\).

The subsequent background has a stiff density and a radiation
density. The radiation fraction is \(f=[1+\exp(-\sigma(N-N_R))]^{-1}\),
where \(N_R=25.0037\) is the original nominal radiation match. The
effective \(\epsilon_H=2.4-0.4f\) gives
\[
 \ln(H/H_R)=-2.4(N-N_R)
       +\frac{0.4}{\sigma}\ln[1+\exp(\sigma(N-N_R))].
\]
This preserves the stiff and radiation asymptotes. \(H_R\) is the
radiation asymptote evaluated at 200 GeV and fixes the same comoving
entropy. The actual radiation temperature grows and then falls
during continuous conversion, with maxima of 179 and 172 GeV in the
two histories. The transfer \(Q\) in
\(\dot\rho_s+4.8H\rho_s=-Q\) and
\(\dot\rho_r+4H\rho_r=Q\) is
\[
 Q=H\rho_{\rm tot}(\sigma-0.8) f(1-f)>0 .
\]
The sum conserves total stress. Direct finite differences of the
background obey \(d\ln H/dN=-\epsilon_H\) through both seams
to \(5.2\times10^{-6}\) absolute. A rolling scalar realization of the
stiff fluid during conversion has \(\dot\chi^2=1.6\rho_s\),
\(V(\chi)=0.2\rho_s\), and the positive effective decay coefficient
\(\Gamma/H=(\sigma-0.8)f/1.6\) along its trajectory. This is an effective
energy-exchange prescription with adiabatic perturbations. It
supplies no decay Lagrangian or new direct vector interaction.

| Effective history | Narrow | Wide |
|---|---:|---:|
| Inflation end blend half width \(\delta\) in e-folds | 0.06 | 0.15 |
| Radiation conversion rate \(\sigma\) per e-fold | 24 | 12 |
| Physical pivot-to-end e-folds | 60.6424 | 60.6252 |
| Integrated vector rest density / frozen target | 0.999972 | 0.999966 |
| Vector energy or pressure fraction through reheating, upper estimate | \(2.45\times10^{-4}\) | \(2.30\times10^{-4}\) |
| Tensor energy or pressure fraction through early radiation, upper estimate | \(4.89\times10^{-4}\) | \(4.59\times10^{-4}\) |
| Extra radiation \(\Delta N_{\rm eff}\), upper estimate | 0.001092 | 0.001026 |
| Uncomputed hard energy / computed hard energy, upper bound | 1.49% | 0.00681% |

The pivot exits long before either blend at the same physical scale
and Hubble rate. The physical end is where the smoothed trajectory
reaches \(\epsilon_H=1\), which changes the e-fold count slightly.
The upper bound on the vector stress through the end of the early
radiation spectator window is \(0.00863\), including the increasing
cold rest density. The audit ends the spectator window at a 1%
standard matter fraction. The later vector matter is the intended
dark matter. The rest-density entries include modes through \(q=10^8\).
The positive tail can add less than \(4.1\times10^{-4}\) of the target
under the inflationary spectrum and the hard-mode controls below.

## How the stress and spectrum were checked

For the soft longitudinal modes use the full Proca energy and
isotropically averaged pressure. With \(K=a^2m^2/(k^2+a^2m^2)\)
and \(P=K A_L'\) in physical conformal units, their quadratic
numerators, with a common factor \(1/(2a^4)\), are
\[
 E_L=|P|^2/K+a^2m^2|A_L|^2,\qquad
 p_L=\big[(3-2K)|P|^2/K-a^2m^2|A_L|^2\big]/3 .
\]
The mode equations give \(E_L'-4\mathcal H E_L+
3\mathcal H(E_L+p_L)=0\). Retain the amplified growing component
and use the frozen tensor spectrum to normalize the soft vector
distribution. The pressure can be negative while a mode is frozen.
Use the positive stress envelope for the gap and the unresolved
momentum bands. Transverse Proca modes have a conformal massless
limit. Their finite-mass gravitational production is subleading
in this hierarchy.

For \(K_h=k/(a_{\rm end}^{\rm sharp}H_{\rm end}^{\rm sharp})\) from
0.001 to 64, evolve the vacuum modes through the inflation-ending
blend and project onto the exact stiff positive-frequency
Hankel basis. Reheating follows with \(k/(aH)>10^{12}\) throughout
the hard band, so its occupation transfers adiabatically.
The Bessel continuation below this band agrees with the exact
growing-mode normalization to \(3.6\times10^{-6}\) on the overlap.
The early time traces use fourth-order adiabatic UV subtraction
with an IR window. The local curvature terms are kept separate
under the spectator GR renormalization convention. Their scale is
below \(1.6\times10^{-9}\) of the background through the blend.
The quoted short-wave tensor stress counts both helicities.
Outside the horizon the quadratic mode stress is an envelope
estimate, with the inflationary contribution of order \(H^4\).
Late pressure averages over rapid phases within a Hubble time.
The archived time traces and envelopes cover the inflation-ending
blend, the stiff evolution and reheating. State-dependent stress
stays small enough to leave the fixed background consistent at
the registered tolerance.

The hard energy integrand is \(K_h^4|\beta_{K_h}|^2\).
The narrow history places 90% of its resolved hard energy between
\(K_h=0.011\) and 29, and the wide one between 0.010 and 12.
About 41.3% and 38.0% of the resolved hard energy lie above
\(K_h=1\). This is sensitive to the inflation-ending transition.
The rest density instead receives 46.4% from
\(0.25\le q=k/k_\ast\le 4\), far below those hard scales.
Tightening the hard grid and vacuum start changes the hard
integrals by less than \(4.9\times10^{-7}\) relative. This
numerical check is separate from the transition sensitivity and
the UV bound.

The first narrow UV estimate using numerical derivatives missed
the 5% gate at 7.32%. It is retained in the audit. Differentiating
the underlying background equations and blend analytically controls
the fourth-order residual. For an adiabatic frequency
\(W=k-U/(2k)+(U''-U^2)/(8k^3)\), the residual variation of constants
gives \(|\beta_k|\le\sinh(C/k^5)\) on the hard tail. The bound
integrates the absolute residual with a factor-two quadrature
margin and uses two quadratures of the derivatives. An independent
oscillatory estimate from the second-order residual also meets the
energy gate. These bounds and the analytic low-\(K_h\) continuation
avoid a cutoff at the inflation horizon. They apply to the stated
smooth effective trajectories.

## Radiation and cold matter

At 1 MeV the horizon is \(k_{\rm BBN}/k_\ast=6.36\times10^{-7}\).
The hard tensor and vector energy are well inside it. Include the
growing-mode contribution down to that horizon and count even the
hard vector energy that will later become nonrelativistic as
extra radiation. With the specified thermal entropy, a relativistic
fraction \(R\) at early radiation converts to
\(\Delta N_{\rm eff}=2.8588 R\).
In the narrow history the bounds are
\(R_{\rm GW}\le2.55\times10^{-4}\) and
\(R_{\rm vector,rel}\le1.28\times10^{-4}\).
The resulting \(0.001092\) is below the registered 0.1 gate
and far below the precision of the BBN-only abundance analysis
\(\Delta N_{\rm eff}=-0.10\pm0.21\) using PDG abundances and
PArthENoPE. The conservative radiation estimate is also small
compared with the Planck+BAO \(N_{\rm eff}=2.99\pm0.17\) constraint.
These are radiation checks. The prior BK18 and Planck tensor/tilt
marginals were already inspected before this audit and cannot be
reserved observations. No new tensor or tilt holdout was used.
[BBN Table 2](https://arxiv.org/abs/2401.15054).
[Planck+BAO](https://arxiv.org/abs/1807.06209).

Using \(v=q/\sqrt{x^2+q^2}\), \(x=a/a_\ast\), and integrating
a conservative radiation conformal time through equality, the
0.01 Mpc streaming cut admits \(q\le8.36\times10^6\).
Those modes also have \(v_{\rm eq}<3.4\times10^{-6}\).
The positive tail and the computed hard modes bound the rest
density outside both cuts to 0.116%. The cold condition thus
holds for at least 99.88% of the final rest density. This is a
kinematics test with the prescribed spectrum. A small-scale
structure prediction would require the density correlations
and nonlinear evolution.

## Reproduce

With Python, NumPy, SciPy and SymPy in this directory, run

```sh
OPENBLAS_NUM_THREADS=1 python stiff_consistency.py narrow
OPENBLAS_NUM_THREADS=1 python stiff_consistency.py narrow --fine
OPENBLAS_NUM_THREADS=1 python stiff_consistency.py wide
OPENBLAS_NUM_THREADS=1 python stiff_consistency.py wide --fine
python stiff_ledger.py
```

[Audit result](results/stiff_consistency.json).
[Initial UV bound](results/consistency_uv_attempt.json).
[Hard mode and vector transfer code](stiff_consistency.py).
[Envelopes and radiation calculation](stiff_ledger.py).

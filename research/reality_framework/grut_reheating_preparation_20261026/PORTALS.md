# Portal conventions, sources, and cold-core matching

This file derives source ledgers, not a late hidden-plasma relic. The
original e_X=1e-7 prompt-radiation heavy point supplies the primordial
vector band. Its parameters are pinned in SOURCES.json. None of that
band is reinterpreted as a computed vector population in the new inflaton
histories.

## Kinetic mixing and widths

Use `-epsilon_BX B_mn F_X^mn/2`. Above electroweak breaking the induced
vector coupling is `epsilon_BX g_Y Y`, to leading order in small mixing.
The current fields have no simultaneous hypercharge and hidden charge,
so the one-loop mixed charge trace vanishes. This does not set a freely
supplied matching coefficient to zero. An exact hidden charge-conjugation
symmetry could forbid it, provided the complete holder/CP sector also
respects that symmetry.

With massless SM final states, sum of Weyl-fermion `Y^2=10` and sum of
complex Higgs-scalar `Y^2=1/2` give

\[
\Gamma_{\rm SM}=\frac{41\epsilon_{BX}^2g_Y^2m_A}{96\pi},\qquad
\Gamma_X=\frac{e_X^2m_A}{48\pi}
(1-4m_X^2/m_A^2)^{3/2}.
\]

Declare `g_Y=.37`, `g*=g*s=106.75` and constant
`m_A=3.2003908656e6 GeV`, so `Gamma_X=2.1223251087e-10 GeV`.
Define `B=Gamma_X/(Gamma_X+Gamma_SM)`. A cold primordial vector population
in radiation domination decays earlier by the total width; its hidden
daughter conformal energy is therefore `R_0 B^(3/2)`, with
`R_0=40.5626228853`. This uses the already calculated finite positive band,
not all possible initially independent quantum modes.

## Thermal inverse decays cannot be omitted

Let `x=m_A/T`, `H=H_m/x^2`, `H_m=sqrt(pi^2 g*/90)m_A^2/M_r`,
`C=45g_A/(4pi^4 g*s)` and `g_A=3`. Then
`Y_A,eq=C x^2 K_2(x)` and the MB averaged inverse Lorentz factor is
`k=K_1/K_2`. The numerical moment closure evolves

\[
\begin{aligned}
\frac{dY_A}{d\ln x}&=\frac{x^2 k}{H_m}
[\Gamma_{\rm SM}Y_{A,eq}-(\Gamma_X+\Gamma_{\rm SM})Y_A],\\
\frac{dY_{\rm pair}}{d\ln x}&=\frac{\Gamma_X}{H_m}x^2 kY_A,\\
\frac{dJ}{d\ln x}&=\frac{\Gamma_X}{H_m}x^3Y_A,
\qquad J=\frac{\rho_{X,\rm relativistic}}{sT}.
\end{aligned}
\]

The energy equation uses the exact vacuum decay identity
`Q_energy=Gamma_X m_A n_A`. The number-time profile uses an assumed MB
kinetic shape at the bath temperature, not a separately proven
thermalization of A. The late nonrelativistic limit is controlled; a full
momentum-dependent vector calculation is absent. Start at x=1e-4 with no
thermal A, end after depletion and integrate both moments. No negative
population is clipped: negligible negative numerical tails are monitored
against the solver's scaled absolute tolerance.

There is a useful exact integrated identity within this linear source
problem, independent of the intermediate time dilation:

\[
Y_{\rm pair}(\infty)=C\frac{3\pi}{2}
\frac{\Gamma_X}{H_m}(1-B),
\quad \int_0^\infty x^3K_1(x)dx=3\pi/2.
\]

Every produced on-shell A eventually decays, with probability B into a
pair. Its comoving energy cannot decrease while it freely redshifts with
fixed nonzero mass: `a E=sqrt(q^2+a^2m_A^2)`. Immediate decay therefore
gives a lower bound on daughter conformal energy. Using
`integral x^4 K_2 dx=15pi/2`, this is

\[
J\ge5Y_{\rm pair}.
\]

It does not require A to maintain the bath's kinetic shape between
production and decay. It does require the stated fixed-mass linear system,
constant relativistic degrees of freedom, no uncounted scattering energy
loss or significant hidden inverse reaction, and subsequent daughter
energy counted before hidden processing.

For the implemented closure the positive delay can also be written
`B integral Y_A d[x K_2/K_1]`. The weight increases with x: `x K_2/K_1`
is the mean energy in units of T under the inverse-decay energy density
proportional to `sqrt(epsilon^2-x^2)exp(-epsilon)`, epsilon>=x. Increasing x
raises this mean by monotone likelihood ratio. The code checks monotonicity
on 2000 log-spaced nodes as well as the two analytic Bessel integrals.

Normalize conformal energy with the old
`E_X=38386.250633... H_i^4` at a=1 and define
`F_R=4(M_r/H_i)^2/E_X`. Thus `R_thermal=F_R J`. We find

\[
R_{\rm sat}=5C\frac{3\pi}{2}\frac{\Gamma_X}{H_m}F_R
=39183.973956,
\qquad
R_{\rm total}\ge R_0 B^{3/2}+R_{\rm sat}(1-B).
\]

The right hand side is minimized at B=1 because
`R_sat > 3R_0/2`. Its minimum is R_0. At epsilon_BX=1e-7 the computed
time-delay moment gives `R_total=983809.9`, and
`Y_pair=1.66294e-7`. At large mixing the source approaches the prompt
bound. The 19-point scan, tighter tolerance, analytic number and quadrature
controls are diagnostic checks, not a viable-point scan.

Thermal masses/statistics/widths, electroweak crossover, off-shell SM
scatterings, changes to the vector mass, inverse X reactions, hidden
collisions and the final f_sym are not solved. The source's maximum energy
is small relative to the hot SM radiation in these runs; it can nevertheless
be large relative to the very dilute inherited hidden condensate.
An altered reheating or vector spectrum needs new matching. An abundance
solver must retain source time and momentum information, not replace the
hard daughters by a thermal population.

## Higgs portal and matching

Use `V contains b phi^2 H^dagger H+lambda_HX |X|^2 H^dagger H`, with
vacuum subtractions/mass matching as in the older portal packet.
The earlier `kappa_Hphi=2b`. For `H=(0,(v_H+h)/sqrt2)` and small mixing,

\[
\theta\simeq\frac{2b v v_H}{m_h^2-m_s^2},\quad
\Gamma_{h\to ss}=\frac{b^2v_H^2}{8\pi m_h}\sqrt{1-4m_s^2/m_h^2},\quad
\Gamma_{h\to X\bar X}=\frac{\lambda_{HX}^2v_H^2}{16\pi m_h}
\sqrt{1-4m_X^2/m_h^2}.
\]

These displayed widths are leading contact-channel estimates; additional
mixing amplitudes and finite-temperature effects are not a complete decay
calculation. At mh125GeV, v_H246GeV and the declared SM width4.1MeV,
the archived ATLAS .107 invisible fraction gives b_max=.0050501 and
lambda_HX,max=.0071419, with one new channel at a time. Production is
assumed SM-like and the final states invisible. No claim of the latest or
complete bound is made.

In a vacuum-mass, constant-g*, thermal-h control, `b=5e-11` gives
`Gamma_hss=4.81568e-20 GeV`. If each s promptly yields an X pair, then

\[
Y_{\rm pair}=2\frac{135}{8\pi^3g_{*s}}
\frac{\Gamma_{hss}}{H(T=m_h)}=2.23656\times10^{-8}.
\]

This is an unprocessed source, not a predicted relic. Thermal phi/X
masses, phase restoration or subsequent reactions can change that result.

For the isolated H-X bridge, the two-real-component X trace is
`Tr M_X^4=2(m_0^2+g_X^2 phi^2+lambda_HX H^dagger H)^2`.
The cross term and scale invariance of the one-loop potential give

\[
\left.\beta_b\right|_{g_X^2\lambda_{HX}}
=\frac{g_X^2\lambda_{HX}}{8\pi^2}.
\]

This is one beta contribution at scales where these modes are active,
not a full RG system or finite threshold matching.

## Cold soliton leakage and scattering

Integrating out the two real Y components with the broken-phase mass
splitting gives, in the V convention,

\[
V_{\rm eff}\supset\alpha_6 X^6+\mathrm{h.c.},\quad
\alpha_6=\frac{\kappa^2\mu h_\Sigma}
{\sqrt2(M_Y^4-2\mu^2h_\Sigma^2)}.
\]

At r=sqrt2 mu h_Sigma/M_Y^2=.45, the exact/leading coefficient ratio is
`1/(1-r^2)=1.2539185`. The old included cold coherent and vacuum-Born
antiparticle-pair rates scale as `|alpha_6|^2`, giving `Gamma t_0=7.08368e-17`.
It is inappropriate to replace this by a free-gas n_X^6 reaction law for
the coherent bound background. Additional channels and the full gauged
in-medium fluctuation problem remain open.

For a fixed low-energy heavy-photon contact, the zero-q point-source
benchmark is

\[
\sigma_A(0)=\frac{\mu_A^2}{\pi}
\left(\frac{Q Z e_X e_{\rm EM}\epsilon_\gamma}{m_A^2}\right)^2.
\]

The three target estimates use the declared fixed flux, all coherent form
factors set to1, and an epsilon_gamma=1 reference. Physical small-mixing
rates scale as epsilon_gamma^2; this is not the validity claim of a
large-mixing diagonalization. The inherited tree soliton-Higgs identity
gives the scalar benchmark in REPORT.md; no finite-q data fit is performed.

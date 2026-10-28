# Canonical vertices, three channels, and what capture means

28 October 2026. Keep all GRUT-M1 inputs fixed. In this packet m_X=1,
z=sqrt(2)X/v and S=phi/v; radial energy and charge from the inherited engine
are multiplied by g_X² relative to physical E/m_X and Q. The positive core is
z=F(r) exp(-i Omega t), with Q=200000, Omega=.88172049195 and E/(Qm)=.93739999053.
These last two numbers are different quantities: dE/dQ=Omega*m_X, not E/Q.

## Vacuum normalization

With h=phi-v, the mixed potential is

\[
g_X^2(v+h)^2|X|^2.
\]

The canonical hXXbar vertex is 2g_X²v, and hhXXbar is 2g_X². The neutral vertices
are g3=6(lambda_phi-epsilon)v and g4=6lambda_phi. Hence, for distinguishable
X and Xbar daughters,

\[
\Gamma_{h\to X\bar X}
=\frac{(2g_X^2v)^2}{16\pi m_\phi}
\sqrt{1-\frac{4m_X^2}{m_\phi^2}}
=\frac{g_X^4v^2}{4\pi m_\phi}\sqrt{\frac13}.
\]

The value is1.73327704e-5 GeV=1.87565899e-4 m_X; lifetime3.7975e-20 s, or
5331.46/m_X. It is20 times the dispatch's numerical9.38e-6m_X. The dispatch's
first displayed formula is already low by2; its following substitution loses
another factor10 at g_X=.1. Sum of the two identical-real-component decay widths
and direct angular two-body phase space give the same answer. The earlier
e71730c kinetic input already has this exact correct width and vertex. M1's
changed neutral cubic does not change this leading decay width.

The standard two-body normalization is in the [PDG kinematics review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf),
Eqs49.11 and49.18. This is a canonical-action calculation, not a new measured width.

## A sub-vacuum mass term does not determine capture

The core's bare S*m_X term alone is not its fluctuation mass. Define
a=lambda_X/g_X²=.1, b=lambda_phi/g_X²=3, c=epsilon/g_X²=.2. Linearize

\[
z=e^{-i\Omega t}\{F+u(r)e^{-i\nu t}+v^*(r)e^{i\nu t}\},\qquad
\delta S=w(r)e^{-i\nu t}+w^*(r)e^{i\nu t}.
\]

For y=r(u,v,sqrt(2)w), the three radial equations are

\[
-y''+\left[
\begin{pmatrix}
S^2+2aF^2&aF^2&\sqrt2FS\\
aF^2&S^2+2aF^2&\sqrt2FS\\
\sqrt2FS&\sqrt2FS&D
\end{pmatrix}
+\frac{\ell(\ell+1)}{r^2}I
-\operatorname{diag}\{(\nu+\Omega)^2,(\nu-\Omega)^2,\nu^2\}
\right]y=0,
\]

where D=b(3S²-1)-6c(S-1)+F². The sqrt(2) normalization makes the matrix symmetric.
The quartic supplies both a normal mass correction and an anomalous mode coupling;
it is not a charge-odd static vector force. At the center, S=.12577 while the
normal charged diagonal S²+2aF²=1.40021. The full coupled problem is essential.

The two-frequency rotating-background mechanism has published prior art in
[Saffin, Xie and Zhou, Q-ball Superradiance, PRL131,111601](https://arxiv.org/html/2212.03269v2),
especially Eqs4-9 and14-17. Their main model differs from M1. The real binding
field adds our third channel; the M1-specific profiles and numerical cross
sections here are new calculations, not independent replication of that paper.

## Channel energies and conservation

An incoming antiparticle has E=sqrt(1+k²), so choose nu=E+Omega. The asymptotic
channels are:

| Channel | Energy in m_X units | Charge |
|---|---:|---:|
|Incoming/reflected Xbar|E|-1|
|Outgoing X|E+2Omega|+1|
|Outgoing neutral h|E+Omega|0|

The X channel is always open. The neutral channel requires E+Omega>sqrt(6):
k>1.20743547. A low-energy incoming antiparticle cannot produce a free on-shell
neutral scalar just because the scalar mass term is reduced inside the core.
The in-core field can respond below this threshold, but its neutral exterior
mode is evanescent. Off-shell exchange is not an outgoing particle carrying
energy to infinity.

At leading linear order, flux-normalized outgoing probabilities obey

\[
P_{\bar X}+P_X+P_h=1,\qquad
\Delta Q_{\rm core}=-2P_X-P_h,
\]

\[
E_{\rm out}-E_{\rm in}=\Omega(2P_X+P_h).
\]

This is conservation of E-Omega Q for the perturbations, with the core's charge
and energy change included. In the first process Xbar+Q -> X+(Q-2), the outgoing
X is energetic. At k=0 its energy is2.76344m_X. No heat is silently erased.

After a single free h decay, one antiparticle is returned for every h produced.
Under that explicitly vacuum-decay-only prescription, the surviving antiparticle
fraction is Pbar+Ph=1-PX. Neutral conversion is therefore not automatically an
antiparticle sink. In a plasma the h may scatter, undergo number-changing reactions
or be absorbed before it decays; none of those possibilities is settled by this
vacuum statement. The scattering approximation treats h as a narrow-width channel,
then adds its decay in a separate ledger. Additional quantum/bound-state processes
and resonance damping are not claimed to be complete.

## Emission is powered by excitation

A stable core cannot lower its total ground-state energy by increasing Q and
emitting an antiparticle. For one such emission the cost is approximately

\[
E(Q+1)-E(Q)+m_X\simeq(\Omega+1)m_X>0.
\]

The disturbance supplies that energy. E(Q)/Q can fall while E(Q) itself rises.
The previous hot-core source is therefore not spontaneous ground-state emission
or an energy-free charge pump.

## Numerical normalization

The finite-difference boundary uses the outgoing discrete radial free solution.
For q_i=y_(i+1)/y_i, propagate q_(i-1)=1/[2+h²(ell(ell+1)/r_i²-k²)-q_i]
backwards from a remote free boundary. Im(q_N)/h is the exact open-channel flux
per boundary amplitude squared; closed-channel q is real. One incident v channel
is normalized to unit flux. A real symmetric interior then gives a unitary
three-channel number-flux matrix. Conservation alone is checked against a separate
continuum collocation solution and nonlinear packets; it is not the sole test.

For a plane-wave incident antiparticle, the inelastic partial-wave sums are

\[
\sigma_X(k)=\frac{\pi}{k^2}\sum_{\ell}(2\ell+1)P_X^{(\ell)}(k),\qquad
\sigma_h(k)=\frac{\pi}{k^2}\sum_{\ell}(2\ell+1)P_h^{(\ell)}(k).
\]

An l0 spherical packet measures only its monopole response. It does not by
itself determine the total cross section or an absorption coefficient.

The reverse channel is required too. At the same nu the real symmetric radial
matrix gives a reciprocal flux-normalized S matrix. In particular,
P(Xbar->X)=P(X->Xbar) for a matched angular channel, while the different incoming
momenta imply kbar² sigma_forward=kX² sigma_reverse. The fast outgoing X can
therefore recreate an antiparticle at a subsequent object. The inverse exchange
closes only below incoming X energy1+2Omega=2.76344m_X. Its source X are above
that threshold. At the channel matched to kbar1, sigma_reverse=.30318/m_X²
and the per-particle reverse rate is about.148 times the forward rate at the
same object density. Above neutral opening, X->h followed by h decay is another
antiparticle source. A forward Gamma/H is not the forward-minus-inverse kinetic
sink. INVERSE_RESULTS.json checks full S-matrix unitarity, reciprocity and the
angular cross-section relation directly; no net-erasure history is inferred.

At a weak exterior probe the packet analysis separates directions in Fourier
space. With U=rZ and radial central derivative DZ, use
Dtilde=[(r²-h²)/r]DZ+cos(k_phase*h)Z. Then (U +/- Dtilde/(i k_J))/2 isolates the
two directions, with k_J=sin(k_phase*h)/h and the leapfrog dispersion correction
included. This avoids equating net receiver charge with a reflected particle
number. A no-core packet, small amplitude, phase rotation and grid controls are
reported. The actual packet-frequency distribution is averaged over stationary
responses before comparison.

## Limits of a rate estimate

For an explicitly assumed cold Q200000 population, nobj=eta Y_Q s_SM/Q and
Gamma=nobj<sigma_X v>. It is a conditional rate, not a derived population or
integrated optical depth. The background, momentum distribution, inverse
reactions, extra sources and changing core charge must evolve together.

The hot-box failure does not prove that no finite-time classical resolution
can converge. Likewise, changing to comoving coordinates does not guarantee a
fast cascade is redshifted away. At the supplied late F=.5 or1 densities,
mapping to M1 Y_Q with illustrative g*=g_s=100 gives H/m_X~3e-14 to7e-14.
Expansion over300/m_X is then only~1e-11. This conditional example is not a
formation-temperature calculation; it demonstrates why a physical clock is
needed before using Hubble friction as a numerical cure. A kinetic split must
transfer energy and charge explicitly and avoid counting the same process twice.
In a fixed comoving grid both physical k and its physical Nyquist cutoff scale
as1/a; redshift alone does not change their ratio. Expansion may change the
interaction/cascade rate, but this must be demonstrated with the physical clock.

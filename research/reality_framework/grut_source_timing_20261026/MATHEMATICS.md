# Equations and scope

The homogeneous19-variable Jordan-frame equations, preparation, scalar Hessian and Coulomb-gauge elimination are inherited unchanged frome160356 and2cbc7f1. The new C implementation propagates one12-component preparation-phase tangent vector per momentum instead of all12 columns of a transfer matrix. It is explicitly adapted from modes.c; it is not an independent derivation of the same solver.

For canonical perturbations q,p, the inherited Hamiltonian is

\[
H_2={p^2\over2}+{1\over2}q^T(k^2+M)q
+{e^2\over2k^2}(g\cdot p-g'\cdot q)^2,
\quad
M=a^2V_{rr}+a^2\xi_XR P_X-{a''\over a}I.
\]

Here chi=a r, g is its charge-weighted gauge tangent, prime denotes dimensionless conformal time, and d eta/d tau=1/a. The six-real-field mode includes the longitudinal degree of freedom; metric, inflaton and radiation perturbations are absent. The neutral binding field is not backreacted in this extension. The original code's Ward projection is retained to avoid subtracting large nearly equal gauge terms.

For the same correlated phase preparation as the old packet,

\[
v_\theta(0)=(\partial_\theta\chi,\partial_\theta\chi),
\quad {\delta |X|\over |X|}=\delta\theta\,
{\chi_X\cdot\delta\chi_X\over |\chi_X|^2}.
\]

Physical initial velocities vanish, a=H/H_i=1, and the seed current vanishes, so these initial conformal canonical momenta are consistent with A0 elimination. A finite difference of the same preparation at phase offsets±1e-5 supplies the derivative. Local charge response uses the canonical current q1 p2−q2 p1; the old phase archive uses its negative, and the verification converts that sign explicitly. No sign is selected by fitting the output.

Each momentum is a separate cosine perturbation. Maxima across momenta are not a stochastic RMS or one simultaneous perturbation field. When a seeded linear radius modulation reaches1, that linear solution has certainly lost small-perturbation control by that sampled time. It does not give the true nonlinear breakup time or even guarantee the first continuous crossing was detected: narrow excursions between saved samples can precede it. Responses after each seed's own threshold are mathematical extrapolations only. The .1 threshold is an earlier warning, not a universal failure criterion.

## Decay clocks

Use tau=H_i t measured from the supplied a=1 start. A parent available then, at rest, with fixed late width Gamma has

\[
N(\tau)=N_0 e^{-\gamma\tau},\quad\gamma={\Gamma\over H_i},
\quad F_N(\tau)=1-e^{-\gamma\tau}.
\]

Thus tau_f=−log(1−f)/gamma. This fraction is independent of N0. Relativistic time dilation lowers the early number of decays; the rest-clock value is an upper bound there. The late benchmark decays are nonrelativistic. Finite mass variation, in-medium widths and changing parent populations are not newly computed.

The quiet reheating histories approach

\[
C=a^2H/H_i,\quad a^2=2C(\tau-\tau_*).
\]

C is measured, not set to1. At tau_m=800 we match tau_*=tau_m−a_m²/(2C). The actual numerical history is used earlier. The parent decay clocks can therefore be compared without transplanting the old instant-radiation populations or R normalization into the new history.

For daughters relativistic at injection, the total comoving energy per decay is proportional to a. The prescribed rest-parent source's cumulative conformal-energy clock is

\[
F_E(\tau)={\int_0^\tau\gamma a(s)e^{-\gamma s}\,ds
\over\int_0^\infty\gamma a(s)e^{-\gamma s}\,ds}.
\]

This is injected source energy under the benchmark law, not thermalized heat or energy remaining at a later matching time. For the matched radiation tail,

\[
\int_t^\infty\gamma\sqrt{2C(s-\tau_*)}e^{-\gamma s}ds
=\sqrt{{2C\over\gamma}}e^{-\gamma\tau_*}
\Gamma\!\left({3\over2},\gamma(t-\tau_*)\right).
\]

Early integration uses Gauss quadrature on each PCHIP interval of a(tau); orders4 and8 agree to3.4e-16 in the checked fraction. An independent improper quadrature checks the tail formula. In the late-decay limit, mean a≈sqrt(pi C/(2gamma)); the earlier source convention had C=1. The normalization N0 and the source energy R must still be recomputed on a selected reheating/inhomogeneous history before abundance inference.

## Auxiliary Floquet spectrum and wavelength diagnostic

The earlier isolated complex quartic orbit is closed after one radial period by undoing its O(2) precession. New infrared and band-edge samples retain that same orbit. They confirm mu/k≈.440 as k→0, a primary sampled upper edge3.135–3.140, and a second sampled band5.57–5.98. These are properties of the auxiliary orbit, not constant Floquet exponents for the nonperiodic coupled history. Exactly k=0 remains a symmetry/zero-growth limit.

For a homogeneous comoving charge N_X=a³ n_X/H_i³, the charge in a cube with side one physical wavelength is exactly

\[
Q_{\rm cube}=|N_X|\left({2\pi\over k/H_i}\right)^3.
\]

It is independent of a. At the primary peak this is about1.34–2.09e5 on the two histories, with a much smaller cube for the second band. This is a geometric diagnostic, not the charge of an actual fragment. Multiple unstable wavelengths, their initial covariance, nonlinear transfer, merging and cooling remain uncomputed. No object charge is inserted from this formula.

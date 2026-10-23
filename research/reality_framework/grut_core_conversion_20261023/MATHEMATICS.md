# Charge conversion, visible emission and conserving transport

The action and C1-HP-03/v1 inputs are unchanged. Natural units are used.
The canonical GRUT-RAI repository is not modified or represented as having
certified this calculation. These are author derivations and computations.

## 1. An allowed hidden outlet

The specified U(1)_X is **global**. There is no X gauge field or massless
hidden photon in this action, hence no associated Coulomb focusing.
Gravitational and scalar forces are different interactions.

For a stable ground object, the previous inequality correctly forbids
object(Q)+anti-X -> object(Q-1)+one free on-shell s at the stated cold
energy. It does not exclude

\[
\boxed{\mathrm{object}(Q)+\bar X(E)
\longrightarrow\mathrm{object}(Q-2)+X(E'),\quad
E'=E+E_Q-E_{Q-2}\simeq E+2\omega_Q.}
\]

Net charge is Q-1 on both sides. For the test object,
omega_Q=Omega*m_X, Omega=0.8611951308217225. A cold incoming particle
therefore has E' approximately 251.573 MeV. The outgoing field carries
the incident energy and the removal energy of two core charges. A free
scalar need not be produced. Finite-Q recoil and the change of the profile
are omitted in the static-background rate; this is the leading large-Q,
weak-incident-wave approximation, not an exact single-quantum S matrix.

This is the established Q-ball charge-conversion / energy-extraction
mechanism, specialized here to the displayed quartic two-field action.
It is not a new GRUT principle. The derivation below includes the X
quartic which is absent from some minimal FLS examples.

## 2. The actual coupled equations

Let tau=m_X t, rho=m_X r, X=v psi/sqrt(2), phi=v S, and
psi_0=F(rho) exp(-i Omega tau). With a=lambda_X/g_X^2=.1 and
b=lambda_phi/g_X^2=3, the equations are

\[
\ddot\psi-\nabla^2\psi+(S^2+a|\psi|^2)\psi=0,
\qquad
\ddot S-\nabla^2S+[b(S^2-1)+|\psi|^2]S=0.
\]

Take a real, angular-normalized spherical harmonic and write

\[
\delta\psi=e^{-i\Omega\tau}
 [u(\rho)e^{-i\nu\tau}+v^*(\rho)e^{i\nu\tau}]Y_{\ell m},
\quad
\delta S=[w(\rho)e^{-i\nu\tau}+w^*(\rho)e^{i\nu\tau}]Y_{\ell m}.
\]

The charged linear variation includes
(S^2+2a F^2) delta-psi + a psi_0^2 delta-psi* + 2S psi_0 delta-S.
The real variation includes [b(3S^2-1)+F^2] delta-S and
2S Re(psi_0* delta-psi). Define the radial vector
U=rho*(u,v,sqrt(2)w). It satisfies

\[
U''+\left[K^2-\frac{\ell(\ell+1)}{\rho^2}I-V(\rho)\right]U=0,
\quad
K^2=\operatorname{diag}((\nu+\Omega)^2-1,
 (\nu-\Omega)^2-1,\nu^2-6),
\]
\[
\boxed{V=\begin{pmatrix}
S^2-1+2aF^2&aF^2&\sqrt2 SF\\
aF^2&S^2-1+2aF^2&\sqrt2 SF\\
\sqrt2 SF&\sqrt2 SF&b(3S^2-1)+F^2-6
\end{pmatrix}.}
\]

The potential is real symmetric. The radial Wronskian
Im(U^dagger U') is constant; regularity makes its value zero at the origin.
Open channels use unit-flux incoming and outgoing Riccati-Hankel waves.
Closed channels use the decaying modified spherical Bessel logarithmic
derivative, retaining the virtual real-scalar response.

The numerical implementation solves both a regular-origin radial system
and, for a check, the alternative U=rho^(ell+1) f system with its singular
derivative handled analytically. The radial boundary at rho_min=1e-4 uses
U'=[(ell+1)/rho_min I+rho_min M(0)/(2ell+3)]U, the regular series through
the displayed order. Domain, mesh, tolerance, free-wave and reciprocity
controls are stored rather than inferred from a solver-success flag.

For incident anti-X with E/m_X=e, set nu=e+Omega. In the source band
1<=e<=sqrt(1.5), both charged channels are open and s is closed.
Let P_+ and P_- be outgoing particle-number flux probabilities. Then

\[
P_++P_-=1,\qquad
\Delta Q_{core}=-2P_+,\qquad
E_{out}-E_{in}=2\omega_Q P_+.
\]

Above the scalar threshold there is a third probability P_s and
Delta Q_core=-2P_+-P_s. Flux conservation still holds with all three
channels. One explicitly three-open-channel control tests that case.

Partial-wave cross sections and Moller relative rates against a heavy
stationary object are

\[
\sigma_{-\to+}^{(\ell)}
=\frac{\pi(2\ell+1)}{p_-^2}|S_{+-}^{(\ell)}|^2,
\qquad v_-=p_-/E_-,
\]
\[
p_-^2\sigma_{-\to+}=p_+^2\sigma_{+\to-},\qquad
\frac{(\sigma v)_{+\to-}}{(\sigma v)_{-\to+}}
=\frac{p_- E_-}{p_+ E_+}.
\]

The rate tends to a finite nonzero value at low incoming momentum in the
worked s wave. Higher partial waves and resonances materially affect the
newborn warm source; geometric cross section alone is inadequate.

## 3. A visible channel can be computed without inventing a storage time

At leading portal order an incoming wave can convert into a timelike
Standard Model pair through a virtual s/Higgs response. The Higgs mixing
insertion is spatially proportional to S(rho), so replacing the whole
object by a free-s vacuum width is insufficient. Retain that form factor.

For an incident partial-wave Hankel coefficient B_in, write
Z=rho*sqrt(2)w/B_in. For dimensionless spatial momentum q define

\[
\mathcal F_\ell(q)=\int_0^\infty d\rho\,\rho\,
j_\ell(q\rho)S(\rho)Z(\rho),\quad
s_q=m_X^2(\nu^2-q^2).
\]

The first-order number-flux loss into SM states is

\[
\boxed{P_{SM}^{(\ell)}=\frac{2}{\pi k_-}
\int_0^\nu dq\,q^2|\mathcal F_\ell(q)|^2
\left(\frac{\kappa v v_H}{m_h^2-s_q}\right)^2
\frac{\operatorname{Im}\Pi_h^{SM}(s_q)}{m_X^2}.}
\]

For a tree-level lepton pair,
Im Pi_h(s)=m_l^2 s/(8 pi v_H^2)*(1-4m_l^2/s)^(3/2)
above threshold, and zero below it. The source band never reaches the
muon or hadron threshold even at zero spatial momentum. A separate
photon-coefficient sensitivity repeats the parent's explicit approximate
coefficient; it is not an exact photon or hadronic matching calculation.

The factor 2/pi follows from the angular-normalized Fourier-Bessel
transform. A direct Parseval control and a separately introduced small
local imaginary-potential fixture check the normalization and optical
loss sign. The fixture is removed from all physical production rates.
The physical integral is first order in kappa^2 and uses the undamped
hidden scattering solution. It is a weak-wave cross section, not the
decay lifetime of an arbitrary trapped core excitation.

## 4. Redshift and inverse conversion belong in the same equation

For the illustrative source, use the parent's imposed scalar-at-rest
daughter energy e0=sqrt(1.5), supplied object charge 200000, original
Y_Q, and an explicitly fixed cold object profile. These are not the
unresolved formation output. Define q=+1 for X and -1 for anti-X and
nu=E/m_X-q Omega. A conversion preserves nu. Expansion gives

\[
\frac{d\nu}{d\log a}=-\left[e_q-\frac1{e_q}\right],
\quad e_q=\nu+q\Omega.
\]

Both channels are available only for nu>=1+Omega. Positive particles
below that threshold cannot convert back into free antiparticles.
At a fixed nu, the collision generator acting on (N_+,N_-) is

\[
\frac{d}{d\log a}\binom{N_+}{N_-}
=\frac{n_{obj}}{H}
\begin{pmatrix}-r_+&r_-\\r_+&-r_-\end{pmatrix}
\binom{N_+}{N_-},\quad r_\pm=(\sigma v)_\pm,
\quad n_{obj}=\frac{s_0Y_Q}{200000\,a^3}.
\]

The exact positive two-state exponential handles rapid forward/inverse
collisions. Conservative remapping handles exact free redshift between
collisions. Mesh/time refinements and a factor-ten rate sensitivity are
retained. The table interpolation is a method approximation, not a
uniformly certified interpolation through every narrow resonance.

The free number is conserved. Net conversion Delta N_+ is accompanied by
Delta Q_core=-2 Delta N_+ and Delta E_core=-2 omega_Q Delta N_+.
Energy plus core chemical energy plus explicit expansion work is checked.
The frozen-profile approximation becomes an uncertainty when several
percent of the core charge is transferred; no exact evolving-population
claim is based on it. D particles, collisions among free particles,
object mergers, nonlinear absorption and a hot formation source are
outside this two-species test.

## 5. What the flux ratio actually bounds

For a nonnegative free outgoing ensemble, E>=m_X(N_++N_-) and
Q=N_+-N_-. If Q>=0, then

\[
0\le N_-\le\frac{E/m_X-Q}{2},\qquad
\frac{N_-}{N_++N_-}\le\frac{1-m_XQ/E}{2}.
\]

This is an upper bound, not an equality. A pure X wave with energy
E_particle>m_X has m_X Q/E=m_X/E_particle<1 and no antiparticles.
Interacting net flux cannot automatically be treated as such an ensemble.
The report separately gives measured net ratios and controlled-band
positive/negative-frequency numbers.

The alternative exterior diagnostic masks the core and absorber, forms
U=rho*psi and P=rho*dot(psi), and sine-transforms each equal-time field.
For omega_k=sqrt(1+k^2), the mode numbers are

\[
N_\pm=4\pi h\left|\frac{\sqrt{\omega_k}\widetilde U
\pm i\widetilde P/\sqrt{\omega_k}}{2}\right|^2.
\]

Their difference exactly recovers the masked canonical charge. Their sum
has a free-field interpretation only under the stated exterior approximation.
Mask variation and the local k^2>10 max|delta m^2| cut are disclosed.
The projected snapshot excludes already absorbed radiation and cannot
be added to a finite-time flux integral as a complete source without a
consistent receiver. It does not silently repair the missing soft band.

## 6. The improved radial integrator

The second-order spatial finite-volume Hamiltonian is first rerun with
fourth-order Hamiltonian time stepping, separating time and space errors.
The independent higher-order spatial scheme uses U=rho*psi and
Z=rho*(S-1), uniform nodal weights 4 pi h and the centered fourth-order
Laplacian. Odd endpoint continuation implements regularity. The discrete
gradient energy is -2 pi h sum U.L4.U, so its force is the exact discrete
Hamiltonian derivative. Each Hamiltonian drift and kick preserves charge.

Velocity-Verlet steps are composed as V(w1 dt)V(w0 dt)V(w1 dt),
w1=1/(2-2^(1/3)), w0=-2^(1/3)w1. The distant numerical absorber remains
a separate symmetric velocity operation with explicit charge/energy
receivers; it is not claimed to be fourth order or a physical dissipator.
Its strength and location are varied directly.

Grid agreement at one finite time is weaker than an asymptotic soliton
or cosmological charge distribution. Richardson extrapolation is accepted
only after the observed error order and boundary controls support it;
nonmonotone grid differences are preserved as a failed extrapolation test.

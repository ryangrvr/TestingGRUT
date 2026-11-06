# X8 · A unitary source and a two-time decay kernel

X8 keeps the semiclassical H55 architecture. The scalar potential
and homogeneous expansion are inherited from X5. A sterile Majorana
field has a mass phase $\theta_I=\phi_{j(I)}/f_I$. The transport
calculation chooses the Standard Model lepton–Higgs Yukawa vertex
with no scalar phase in the original mass-phase frame. In covariant
notation the new matter terms have the form

$$
 \mathcal L_N=\sum_I \tfrac12\overline N_I
    \big(i\not\nabla-M_I e^{i\gamma^5\theta_I}\big)N_I
   -\left(Y_{\alpha I}\overline L_\alpha\widetilde H P_R N_I
      +\text{h.c.}\right).
$$

The calculation uses the inherited diagonal vacuum widths.
A specified flavor matrix, thermal state and state correlations
would be needed for flavor-changing or medium terms. The spatial
control resolves one periodic axis and integrates transverse
momenta with an averaged FRW metric. Its expansion uses X7's
finite-grid metric matching. A general curved-spacetime calculation
requires its own covariant matching. [Axion anomalies](https://arxiv.org/abs/2112.07645)

## The unitary mass operator

Let $P_Q$ project onto the longitudinal modes $-N,\ldots,N$ and set
$T=P_Q\theta P_Q$ and $W=\exp(iT)$. The four-component slab uses the
mass operator and chiral kinetic term

$$
 V_Q=M
 \begin{pmatrix}0&W\\W^\dagger&0\end{pmatrix}\otimes I_{\mathrm{spin}},
 \qquad K_Q=\operatorname{diag}(D,-D),
 \qquad D=(q_x\sigma_z+q_\perp\sigma_x)/a .
$$

$V_Q^2=M^2 I$ at every finite cutoff. In $H_Q^2=K_Q^2+V_Q^2+\{K_Q,V_Q\}$,
the remaining phase-dependent blocks are proportional to $[D,W]$.
The derivative identity
$[D,W]=i\int_0^1 e^{itT}[D,T]e^{i(1-t)T}\,dt$ bounds them by the
spatial derivative of the prescribed phase. This fixes the term
with zero phase derivatives in the large-momentum vacuum expansion.
A spatially uniform phase reduces to the X6 helicity Hamiltonian.

For a variation $\delta T$, diagonalize $T=U\Lambda U^\dagger$ and use

$$
 (\delta W)_{ij}^{U}
  =i e^{i(\lambda_i+\lambda_j)/2}
    \operatorname{sinc}\!\left((\lambda_i-\lambda_j)/2\right)
    (\delta T)_{ij}^{U},\qquad \operatorname{sinc}(z)=\sin z/z .
$$

This is the exact matrix derivative, including coincident
eigenvalues. Its contraction with the same covariance gives the
scalar force. The local rotation integrates the reciprocal scalar
kick along with the orbital evolution. The action checks compare
that force with a variation of fermion energy and the pressure with
a variation of scale at fixed canonical momentum.

## Finite response and kinetic matching

At quadratic order the direct vertex is proportional to $T^2$.
For a cosine of harmonic $r$, the diagonal of $T^2$ counts only
neighbors $n\pm r$ inside the retained band. The old projected
pointwise exponential counted both neighbors in the direct vertex,
while only the inside neighbors entered the virtual term. The new
direct term cancels the unmatched contact contribution exactly
at every frequency in the quadratic kernel.

With $E=E_n$, $E'=E_{n+r}$, $k=2\pi r/L$ and
$s=(k/a)^2-(\omega+i0)^2$, the remaining kernel is

$$
 \Pi_Q(\omega,k)=
   s\,\frac{M^2}{f^2a^3}\sum_{\ell,n\,\mathrm{paired}} W_\ell
   \frac{E+E'}{2EE'[(E+E')^2-(\omega+i0)^2]} .
$$

The independent spectral check uses the direct second derivative of
the finite Hamiltonian and its occupied-to-empty transitions.
For the continuum comparison, one Feynman parameter combines the
propagators. The radial integral to a transverse cutoff $p_\perp$
replaces $D_x^{-1/2}$ by $D_x^{-1/2}-(D_x+p_\perp^2)^{-1/2}$, with

$$
 D_x=M^2+(1-x)(q_n/a)^2+x((q_n+k)/a)^2
                   -x(1-x)(\omega+i0)^2 .
$$

The overall weight is $M^2/(8\pi a L f^2)$. The unbounded periodic
reference sums all $n$ and subtracts the zero-external-momentum
kinetic coefficient. The first slab uses that same matching on all
retained modes. Its paired bubble lacks $r$ of those terms near the
longitudinal boundary, so the remainder contains a $Q^{-1}$ term.

The follow-up uses the average of the two zero-external integrands
on each paired mode,

$$
 Z_r^Q=\frac{M^2}{8f^2 a^3}\sum_{\ell,n\,\mathrm{paired}}W_\ell
            (E_n^{-3}+E_{n+r}^{-3}),\qquad A_r=1-Z_r^Q.
$$

It preserves the homogeneous subtraction at $r=0$. The change in
the action vanishes in the cutoff limit. At high momentum the
symmetrized difference starts with two powers of external momentum
or frequency and gives a $Q^{-2}$ quadratic tail. This result is
for the quadratic kernel. The state preparation and nonlinear
time evolution have separate cutoff errors.

At finite cutoff $A$ is diagonal in scalar Fourier modes. Its
scale dependence is included in the pressure. With
$g=a^{-1}\partial_x\phi$, $v=\dot\phi$, and $A'=dA/d\log a$,

$$
 \rho_\phi=V+\tfrac12\langle v,Av\rangle+\tfrac12\langle g,Ag\rangle,
 \qquad
 P_\phi=-V+\tfrac12\langle v,Av\rangle+
        \tfrac16\langle v,A'v\rangle-
        \tfrac16\langle g,(A+A')g\rangle .
$$

The canonical scalar momentum is $a^3Av$. The longitudinal
translation charge uses its contraction with the scalar gradient
and the fermion momentum from the covariance. The matched fermion
ledger moves the known scalar derivative subtraction to both
fermion stress and fermion work, as in X7 with the Fourier
coefficient in place of the scalar coefficient.
[X7 ledger](../x7_regulator/MODEL.md)

For the startup test, expand the adiabatic projector
$P=P_0+P_1+P_2+P_3+\cdots$. In the instantaneous
occupied-to-empty block $P_{n,+-}=i(\dot P_{n-1})_{+-}/\Delta$.
Idempotency fixes the second-order diagonal blocks to
$X_1X_1^\dagger$ and $-X_1^\dagger X_1$, where $X_1=P_{1,+-}$.
The graph coordinate at third order also includes
$X_1X_1^\dagger X_1$. The code evaluates these derivatives around
the prepared state and applies the extra graph term only above
a chosen gap range. Changing that range tests an initial-state
choice, in addition to the cutoff and timestep checks.

## Constant-rotation decay

For one helicity let $k=hp$, $z=k+\delta$, $\delta=\dot\theta/2$,
and $\Omega=\sqrt{M^2+z^2}$. A positive rotating-frame eigenmode
has chiral weights $(1+bz/\Omega)/2$ for $b=\pm1$.
Undoing the phase rotation gives laboratory frequencies
$E_b=\Omega-b\delta$. The massless two-body phase space gives
an integrated lepton current proportional to the pair momentum.
Its chiral contraction is $E_b-bk=\Omega-bz$. The resulting
partial width is

$$
 \Gamma_b=
  \frac{\Gamma_{\rm rest}}{M}
  \frac{1+bz/\Omega}{2}(E_b-bk)\Theta(E_b-p)
  =\frac{\Gamma_{\rm rest}M}{2\Omega}\Theta(E_b-p).
$$

The two branches give the standard time-dilated width at zero
drive. For helicity opposite the drive one branch closes at
$p>M^2/(4|\delta|)$. The other stays open. The total width
is bounded by $\Gamma_{\rm rest}$ in this stationary vacuum
calculation. The numerical check integrates the lepton current
over the boosted two-body sphere and tests the spin update.

The locally constant replay keeps the full Bloch vector. A jump
damps the parallel occupation about $(M,0,z)/\Omega$ at the total
width and damps transverse coherence at half that width. If a
branch loses occupation $\Delta n_b$, it delivers
$E_b\Delta n_b$ to radiation. The sterile laboratory energy
uses $(M,0,k)$, so the scalar kinetic energy must supply the
difference. At fixed scale, for weight $W_q$, the momentum
impulse completing the collision ledger is

$$
 \Delta P_\phi=W_q\left[\frac{\Delta r_z}{4f}
                      +\sum_b\frac{b\,\Delta n_b}{2f}\right].
$$

The scalar midpoint velocity contracts with this impulse to
cancel the combined fermion and daughter energy change. The same
formula gives the rotating-frame charge ledger, with the signed
branch count tallied for each parent. This is the conservation
completion of the diagonal secular update. The full Standard
Model flavor dynamics need a coupled collision calculation.
The replay measures inherited projection gaps on both resolutions.
Its local-rate indicators compare basis rotation to the gap and
threshold drift and expansion to the bath gap. They diagnose
the local approximation on the tested trajectory.

## Two-time vacuum functional

For a conformally massless empty lepton–Higgs bath, comoving
momentum is $q$, helicity momentum $k=hq$, and the conformal
rotating-frame spinor propagator is $U(\eta,\eta_0)$. The
frequency-resolved source for the archived covariance $C_0$ is

$$
 B_b(\omega)=\int d\eta\,g(\eta)\,
       e^{i\omega\eta+ib\theta(\eta)/2}\,
       e_b^\dagger U(\eta,\eta_0) C_0^{1/2}.
$$

At order $Y^2$ the emitted comoving energy in the smooth exposure is

$$
 E^{\rm bath}_b=\frac{\Gamma_{\rm rest}}{2\pi M}
     \int_q^\infty d\omega\,\omega(\omega-bk)
       \operatorname{tr}\big[B_b(\omega)B_b^\dagger(\omega)\big].
$$

This positive functional retains both coherence components.
For a long stationary exposure it reduces to the on-shell
width above. For a finite window, its frequency resolution
can mix power across the bath threshold even on a frozen
background. X8 runs that frozen control at the same central
state and with the same exposure.

The underlying time-domain vacuum kernel is explicit. With an
exponential convergence factor on frequency above $q$,

$$
 K_b(\tau;\epsilon)=\frac{\Gamma_{\rm rest}}{2\pi M}e^{iq\tau}
     \left[\frac{q-bk}{\epsilon-i\tau}
                 +\frac{1}{(\epsilon-i\tau)^2}\right].
$$

Thus the stationary threshold has a nonlocal time response.
A conserving evolution with a continuously coupled bath must
also handle its dispersive self-energy, the initial
sterile–bath correlations and the energy in the interaction.
The archived sterile covariance and the prescribed exposure
specify the test above. The interacting initial state needs
additional bath information.
[Coherent collision integrals](https://arxiv.org/abs/1006.1929) ·
[Neutrino transport with coherence](https://doi.org/10.1007/JHEP12%282024%29169)

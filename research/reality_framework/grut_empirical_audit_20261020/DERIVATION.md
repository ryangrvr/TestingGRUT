# Independent derivation of the exclusion and distance test

20 October 2026. Author rederivation; no independent researcher or empirical
replication. All inequalities in sections 1–4 concern the displayed leading
classical, minimally coupled action. They are not all-orders quantum bounds.

## 1. Normalization from the action

Let `X=(q1+i q2)/sqrt(2)` with canonically normalized real fields. On a spatial
slice with future unit normal n, define `dot(q)=n^mu partial_mu q`. The current
in Candidate 1 gives

\[
n_Q=-n_\mu j_X^\mu=q_2\dot q_1-q_1\dot q_2.
\]

For `X=f exp(-i omega t)/sqrt(2)`, `n_Q=omega f^2`: there is no extra factor
two. A single X quantum defines unit global charge. Changing that convention
requires changing the numerical yield inversely; it cannot change the energy
constraint.

Write `mX=gX v`, `tau=mX t`, `rho=mX x`, `q1=v u`, `q2=v w`, `phi=v S`,
`a=lambda_phi/gX^2`, `b=lambda_X/gX^2`. Then

\[
Q=\bar Q/g_X^2,\quad E=(v/g_X)\bar E,
\quad n_Q=g_Xv^3\bar n_Q,\quad \rho_h=g_X^2v^4\bar\rho_h.
\]

Thus `E/Q=mX Ebar/Qbar`. The simulated classical Qbar is not a physical
particle count. Natural-unit densities convert to SI by powers of hbar*c.

## 2. Sharp potential coefficient, without a soliton ansatz

For `D=0`, the dimensionless potential is

\[
V=\frac a4(S^2-1)^2+\frac12 S^2 F^2+\frac b4F^4,
\quad F^2=u^2+w^2.
\]

Set `y=F^2>0`, `s=S^2>=0`. The ratio is

\[
\frac{2V}{F^2}=s+\frac{a(s-1)^2}{2y}+\frac{by}{2}.
\]

At fixed s the last two terms have infimum `sqrt(ab)*abs(s-1)`.
Minimizing `s+sqrt(ab)*abs(s-1)` on `s>=0` yields

\[
\mu^2=\min(1,\sqrt{ab}),\qquad
m_* =m_X\min(1,(ab)^{1/4})
=v\min(g_X,(\lambda_\phi\lambda_X)^{1/4}).
\]

For `ab<1`, equality occurs at `S=0,F^2=sqrt(a/b)`. For `ab>1`, the
infimum is the free-vacuum limit `S^2=1,F->0`. At `ab=1` there is a degenerate
minimum. This proof neither assumes nor establishes a finite bound object.
Positive D terms only add nonnegative energy, so D=0 is not necessary for
this particular leading-order inequality.

For either sign epsilon=+1 or -1,

\[
\frac12(\dot q_1^2+\dot q_2^2)
+\frac{m_*^2}{2}(q_1^2+q_2^2)-\epsilon m_*n_Q
=\frac12(\dot q_1-\epsilon m_*q_2)^2
+\frac12(\dot q_2+\epsilon m_*q_1)^2\ge0.
\]

All real-field kinetic energies and spatial gradients are nonnegative on a
spacelike slice. Since `U>=m_*^2 |X|^2`, this proves the pointwise proper-matter
bound `rho_h>=m_* abs(n_Q)`, and consequently

\[
\int\sqrt h\,\rho_h\ge m_*\int\sqrt h\,|n_Q|
\ge m_*\left|\int\sqrt h\,n_Q\right|.
\]

Opposite charges, gradients, moving cores, fission, mergers and evaporation
cannot weaken it at this order. Surface energy makes it stricter. No isolated
profile, numerical component identification or cosmological composite lifetime
enters this proof.

## 3. Expansion, conserved yield and gravitational limits

The same proper-slice inequality holds in curved spatial geometry for this
minimal canonical matter action. In FRW, physical gradients have positive
`A(t)^-2` factors. Expansion does not conserve integrated matter energy, but
the instantaneous inequality still applies. Exact U(1), a closed comoving
domain and homogeneous ensemble averaging give `A^3 n_Q=constant`.
Visible comoving entropy conservation gives `Y_Q=n_Q/s_vis=constant`.

The old preparation equates the primordial hidden yield to the **late** baryon
yield, not to the baryon density before electroweak sphalerons. This equality,
its sign and magnitude are additional initial-state data. For present ratio
`r=abs(Y_Q0)/Y_B0`, baryon mass convention mb and BBN-inferred omega_b,

\[
\omega_{h,\min}=\frac{m_*r}{\bar m_b}\omega_b.
\]

The entropy and CMB-temperature factors cancel in this ratio. Late entropy
injection invalidates a specified *initial*-yield relation unless the changed
present r is explicitly evolved. It does not invalidate the formula at a
separately specified present r. Baryon number must be conserved from BBN onward
to use its inferred present density.

This is a proper matter-energy inequality, **not** a bound on arbitrary ADM
mass, gravitational binding energy, black-hole mass or a quantum expectation
of renormalized stress. Connecting it to a cosmological density requires the
leading Einstein constraint, near-FRW averaging, negligible strong binding,
and no loss of the modeled charge behind an unmodeled strong-gravity region.
Backreaction, nonminimal curvature terms, running/loop changes of U and higher
operators require a new calculation. They are not assigned zero uncertainty
by a small numerical discretization error.

For a spatially flat, expanding FRW background, a nonnegative cosmological
constant and nonnegative other densities,

\[
H_0^2\ge(100\ {\rm km\,s^{-1}\,Mpc^{-1}})^2
\omega_b\left(1+\frac{m_*r}{\bar m_b}\right).
\]

Positive radiation or Lambda strengthen the bound. The hidden pressure need
not be cold dust. Negative Lambda, positive spatial curvature (the negative
curvature contribution to H^2), significant binding/backreaction, altered
gravity, hidden-charge violation or dilution of a fixed initial yield are
outside these assumptions. This does not establish that the postulated
history actually remains in this domain until today.

## 4. A distance envelope that avoids a fitted LambdaCDM expansion history

The published combined maser H0 is inferred using flat LambdaCDM with
Omega_m=.315 taken from Planck (Pesce et al. 2020, Eq. 1). Its geometric disk
distances in Table 1 precede that fit. The earlier description of the H0
comparison as entirely CMB/cosmology independent needs this qualification.

At every redshift in the nearby sample, conserved mean baryon and hidden
charge densities and the same energy inequality imply

\[
H(z)\ge H_{\rm floor}(1+z)^{3/2},\qquad
H_{\rm floor}=100\sqrt{\omega_b(1+m_*r/\bar m_b)}.
\]

For ordinary photon propagation in flat FRW,

\[
D_A(z)=\frac{c}{1+z}\int_0^z\frac{dz'}{H(z')}
\le\frac{2c}{H_{\rm floor}(1+z)}
\left(1-\frac1{\sqrt{1+z}}\right)\equiv D_{A,\max}(z).
\]

This upper envelope does not set an unknown EOS or a free H0 to the observed
answer. Actual gradients, kinetic energy or additional matter only lower D_A
in the stated domain. It uses no dark-matter density or Planck Omega_m prior.
It remains conditional on representativeness of the FRW light path, small
lensing/local-distance corrections, calibrated maser disk inference and a
declared treatment of peculiar velocities.

Table 1 velocities are optical-convention `v_obs=c z_obs` in the CMB frame.
For a line-of-sight peculiar velocity beta, use

\[
1+z_{\rm obs}=(1+z)\sqrt{(1+\beta)/(1-\beta)}.
\]

The protocol states observational-error corners before calculation. They are
sensitivity intervals, not exact posterior support or a joint confidence set.
The envelope increases with z for `0<z<1.25`, encompassing the entire sample.
For a given positive lower distance D, its conditional product constraint is

\[
m_Xr\le\frac{\bar m_b}{\mu}
\left[\frac{1}{\omega_b}
\left(\frac{2c[1-(1+z)^{-1/2}]}{100(1+z)D}\right)^2-1\right].
\]

This is a necessary condition, not a fitted abundance or a sufficient allowed
cosmology. It constrains a product of uncalibrated inputs. It does not identify
the FLS/Q-ball mechanism over another charged massive sector.

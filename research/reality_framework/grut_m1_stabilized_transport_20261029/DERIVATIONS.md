# Definitions and derivations

## Fixed-curvature scalar bound

Use fields in units of `H_i`, `R` in units of `H_i²`, and `V` in units of `H_i⁴`. The heavy spectator potential actually implemented is

\[
V=\lambda_X|X|^4+\lambda_Y|Y|^4+M_Y^2|Y|^2
+\lambda_\Sigma(|\Sigma|^2-h^2/2)^2
+2\Re[\kappa Y^*X^3+u(R)\Sigma^*Y^2]+\xi_XR|X|^2.
\]

Here `u=mu+LOCK R`, `LOCK=(10i−mu)/12`, `M_Y=30`, `h=40`, `lambda_X=.001`, `lambda_Y=lambda_Sigma=1`, `xi_X=−5/6`, with the exact remaining constants read from the pinned M1 point. No flat-space parameters are silently refitted.

Write amplitudes `x,y,s >= 0` and reserve `A=.2 lambda_X`. Minimizing the other X quartic against the kappa term gives

\[
.8\lambda_X x^4-2|\kappa|yx^3\ge
-\frac{27|\kappa|^4}{16(.8\lambda_X)^3}y^4.
\]

Thus `B=lambda_Y−27 kappa⁴/[16(.8 lambda_X)³]=0.5283522987>0` remains. The further square bounds

\[
B y^4-2|u|s y^2\ge-|u|^2s^2/B,
\qquad A x^4-a_x x^2\ge-a_x^2/(4A),
\quad a_x=\max(-\xi_XR,0)
\]

yield, after dropping the positive `M_Y² y²` term,

\[
 V\ge\frac{\lambda_\Sigma h^4}{4}
-\frac{a_x^2}{4A}
-\frac{(\lambda_\Sigma h^2+|u|^2/B)^2}{4\lambda_\Sigma}.
\]

This finite lower bound holds for every fixed finite real R and arbitrary scalar amplitudes/phases. The original quartic form is strictly positive away from the origin after the mixed-quartic inequality; it dominates the cubic and quadratic terms at large total scalar amplitude. A scalar-degree statement is not a statement about the dimension of the curvature operator, the gravitational Hamiltonian, quantum consistency or a UV completion.

For action convention `F R/2 − V_flat − R W`, `F=Mpl²−2W`. Along the trajectory `F/Mpl²=1−2(H_i/Mpl)² dV/dR`; this is checked only along the sampled solution. Arbitrarily large scalar fields can violate `F>0`, outside the regime certified here.

## Why both Option 1B claims need qualification

For the old sextic curvature interaction, its dangerous phase gives a term `−a s x⁶`, where `a>0` at fixed nonzero R. The suggested positive `b x¹²` has the correct four-dimensional mass dimension when `b=xi6² |R|/(2 lambda_Sigma Lambda¹⁰)`. With `b>0`, it can dominate the negative contribution after minimizing over `s`: `lambda_Sigma s⁴−a s x⁶` minimizes as `−const*x⁸`, and `b x¹²` dominates that. At R=0 those curvature terms vanish. This is a bounded fixed-R scalar extension, but it is not the claimed perfect square with the existing Sigma quartic. A square formed from `s²` and `x⁶` has cross term `s²x⁶`, not `s x⁶`. The `|R|` dependence is also nonanalytic at R=0.

The other proposal, `c8 R² |Sigma|² |X|⁶/Lambda¹⁰`, has mass dimension **two**, not four, for dimensionless c8. The denominator would need `Lambda⁸`. Even after that correction, it does not by itself stabilize the old sextic interaction: choose `Y=0`, a fixed Sigma phase maximizing the negative term, and fixed `0<s<a/b`. The X⁶ coefficient `b s²−a s` is negative. Sending X to infinity then beats the positive X⁴ term. This is an explicit unbounded direction, not a numerical failure to find a minimum.

## Reheating and canonical transfer

Real canonical fields are `r=sqrt(2)(Re X,Im X,Re Y,Im Y,Re Sigma,Im Sigma)`. The prescribed dimensionless background is

\[
\dot I=J,\quad \dot J=-(3H+0.5)J-9I,\quad
\dot\rho_r=.5J^2-4H\rho_r,\quad
H^2=(J^2/2+4.5I^2+\rho_r)/3,
\quad R=18I^2-J^2.
\]

Initial `I=sqrt(6)/3`, `J=rho_r=0`, `a=1`, and hidden velocities zero. Scalars obey `rddot+3H rdot+V_r=0`. The hidden energy ledger integrates `dE/dt=−3H rdot²+Rdot V_R`, with `Rdot=54 I J+2(3H+.5)J²`. Charge convention is `N_X=−2a³ Im(X* Xdot)`; total gauge charge weights are `(1,3,6)`.

For conformal variables `chi=a r`, the Coulomb-gauge scalar Hamiltonian gives

\[
q'=p+c g j,\quad p'=-(k^2+M)q+c g'j,\quad
j=g\cdot p-g'\cdot q,
\quad c=e^2/k^2,\quad g=T\chi,
\quad M=a^2(V_{rr}-R/6).
\]

Primes denote conformal time. `T` has blocks `q_charge[[0,−1],[1,0]]`. The Ward identity `M g=a³ T V_r−a² R g/6` is used to avoid subtractive loss in the gauge direction. The `k=0` calculation sets `c=0` and restricts initial scalar data to the five zero-gauge-charge physical eigenmodes; finite `k` tests at `e=1e−6` check that subspace. This is not a calculation of the longitudinal-vector quantum vacuum or of metric/inflaton perturbations.

For a constant de Sitter scalar mass, define `S=sqrt(2k³) delta r/H_i` and `N=ln(aH_i/k)`. Then `S_NN+3S_N+(exp(−2N)+m²/H_i²)S=0`. A positive-frequency adiabatic initialization at `k/(aH_i)=1000` is compared with 2000. For N beyond 8, the negligible-gradient equation is propagated analytically with `exp(−3N/2)` damping and frequency `sqrt(m²/H_i²−9/4)`.

At reheating start the canonical columns are `(U S,U(S+S_N))`, because `a=H/H_i=1`. If L is the late fractional-charge response row, the dimensionless scalar charge power is `sum_a |L column_a|²/(4pi²)`. This fixes the phase-space normalization. The late response is averaged over `640 <= H_i t <= 800`; transfer and direct phase-perturbed trajectories use the same average. The covariance uses five independent scalar modes and an imposed vacuum/duration, not the old unnormalized Gaussian prescription.

## Number-conserving elastic gas evolution

Let `x=a/a_i` with `a_i` at 2 MeV, and `epsilon=x² p²/(2m_X²)`. Redshift is absent from the comoving energy equation. The threshold amplitude is `4lambda_X−g_h²/3=−0.00933333` for `g_h=.2`, and the identical-particle cross section is `sigma*m_X²=amplitude²/(128pi)`.

For incoming comoving energies `epsilon_i,epsilon_j` and incident angle cosine mu, the relative speed in comoving units is `sqrt(2(epsilon_i+epsilon_j)−4sqrt(epsilon_i epsilon_j)mu)`. Isotropic CM scattering makes either daughter energy uniform between

\[
\frac{\epsilon_i+\epsilon_j}{2}
\pm\frac12\sqrt{(\epsilon_i+\epsilon_j)^2-4\epsilon_i\epsilon_j\mu^2}.
\]

Both daughters are deposited with hat functions, exactly preserving number and comoving kinetic energy. Rates sum unordered pairs with the factor one-half on identical grid indices. Incident-angle quadrature and the finite support are separately refined/audited. Unsupported final states are omitted symmetrically and their event probability is recorded.

The exposure satisfies `ds/d ln x=sigma_tilde N_plus (ncore20/m_X³) m_X/(1000 H(x) x⁴)`. The factor 1000 inherits dilution between the old 20-MeV normalization and its `a/a0=10` checkpoint. It is a declared inherited-density assumption, not a re-derived formation density. H and the photon temperature include equilibrium electron/positron entropy transfer and instantaneous neutrino decoupling in the parent background.

The collision count is integrated as `2 events/N`; the viscosity proxy is `2/3` of the isotropic threshold rate. The finite-grid Maxwell stationarity residual is reported, rather than presumed exactly zero. The final fitted phase-space function uses number weighting `q² f(q)dq`; matching its velocity moments does not determine perturbations.

## Conditional density hierarchy

Use conformal Newtonian potentials `Phi=Psi`, `N=ln(a/a_i)`, `K=k/(aH)`, and `u=theta/k`. The gas perturbation is expanded in the usual Legendre hierarchy with moments `Psi_l(q)` and nonrelativistic `v=q exp(−N)` (q here is momentum in units of `m_X` at the initial checkpoint). With `d=d ln f0/d ln q`,

\[
\Psi_{0,N}=-Kv\Psi_1-\Phi_N d+C_0/H,
\quad
\Psi_{1,N}=Kv(\Psi_0-2\Psi_2)/3-K\Phi d/(3v)+C_1/H,
\]

and for `l>=2`, `Psi_l,N=Kv[l Psi_(l−1)−(l+1)Psi_(l+1)]/(2l+1)+C_l/H`.

Number averages use generalized-Laguerre quadrature for `q²f0dq`. Define `delta=<Psi0>`, `ug=<v Psi1>`, `theta_T=<q²Psi0>/<q²>−delta`, and temperature basis `b=−(3+d)/2`, normalized to its discrete energy moment. Our conserving relaxation prescription is

\[
C_0/\nu=-\Psi_0+\delta+b\theta_T,
\quad C_1/\nu=-\Psi_1-d u_g/(3v),
\quad C_{l\ge2}/\nu=-\Psi_l.
\]

It conserves the discrete number, energy and momentum moments; tests use arbitrary nonequilibrium perturbations. For the non-Maxwell fitted shape it is a declared moment-preserving ansatz, not the exact elastic equilibrium operator. The measured viscosity proxy supplies `nu/H`; factor-two changes are controls.

The metric constraints implemented are

\[
\Phi=-\frac{3}{2K^2}\left[\Omega_r\delta_r+\Omega_m\delta_m+
\frac3K((1+w_r)\Omega_r u_r+\Omega_m u_m)\right],
\quad
\Phi_N=-\Phi+\frac{3}{2K}[(1+w_r)\Omega_r u_r+\Omega_m u_m].
\]

Cold matter obeys `delta_c,N=−K uc+3Phi_N`, `uc,N=−uc+KPhi`. Radiation follows the barotropic fluid equations with `w_r` and adiabatic sound speed from the supplied electron/positron background. Analytic radiation-era initial conditions satisfy the Einstein constraints. At formation, gas density and velocity are set equal to cold matter's, `theta_T=(2/3)delta_c` by default (zero in a control), and `l>=2` moments vanish. Already-subhorizon modes inherit this artificial matching too.

Radiation oscillations are averaged after `K=150,300,450` in the respective controls. Quadrature increases from 24 to 48 to 96 nodes; hierarchy truncation from 100 to 220 to 260. The top moment uses zero next multipole. The test gas fraction is the parent dark-matter gas fraction multiplied by `.12/(.12+.0224)`. The background keeps parent `Omega_m=.315`; these rounded background/composition choices and cold baryon treatment are test-model conventions, not a precision cosmological fit. The full formation source, photon/neutrino/baryon hierarchy, exact linearized collision kernel, metric anisotropic stress and nonlinear evolution remain outside this calculation.

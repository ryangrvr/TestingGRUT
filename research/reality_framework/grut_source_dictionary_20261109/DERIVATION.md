# D1: operational source dictionary and its finite-spacing limitation

Historical base: N22 result 748c74ff65206aff27ac550e404b02fab0989f9a. N22 remains a completed finite numerical experiment. This note derives an operational realization of its detector and identifies what that realization does and does not fix. No probabilities enter the derivation.

## 1. Exact electric-energy variable of the shifted lattice Hamiltonian

Let epsilon_i=(-1)^i, N be even, and sum_i Z_i=0. The archived Hamiltonian in physical units is

    H_a = H_hop + (g^2 a/2) sum_{j=0}^{N-2} L_j^2
          + (m_lat/2) sum_{i=0}^{N-1}(1+epsilon_i Z_i),
    H_hop = (1/(2a)) sum_j (T_j+T_j^dagger),
    T_j = sigma_j^+ sigma_{j+1}^-.

Gauss law gives L_j=S_j+1/4+epsilon_j/4, S_j=(1/2)sum_{i<=j}Z_i. Set y_j=L_j-epsilon_j/4. Summing the alternating links gives the exact identity

    sum_j epsilon_j y_j = (1+sum_i epsilon_i Z_i)/4.

Consequently, at m_lat=m_R-g^2 a/8,

    H_a = H_hop + (g^2 a/2)sum_j y_j^2
          + (m_R/2)sum_i(1+epsilon_i Z_i) - g^2 a(N-1)/32.

This is an operator identity within the even neutral sector, including the finite-chain constant. At m_R=0, y is the variable appearing in the smooth electric-energy square. The known mass shift is discussed by [Dempsey et al., Eq. (3)](https://arxiv.org/html/2206.05308v3); the finite-open-chain algebra above is derived directly for this code, rather than inferred from periodic translation symmetry.

## 2. A calibrated background-flux source

Introduce dimensionless link sources q_j into the y-square. The source-dependent difference is

    H_y[q]-H_y[0] = g^2 a sum_j q_j y_j + (g^2 a/2)sum_j q_j^2,
    (g^2 a)^(-1) dH_y/dq_j at q=0 = y_j.

A conventional flux source replacing L_j by L_j+q_j instead produces an additional known scalar g^2 a sum_j epsilon_j q_j/4. Subtracting that source-only term gives the y-source family. This is a physical calibration of the apparatus phase, not an adjustment inferred from observed probabilities.

Use Phi=E/g=phi/sqrt(pi) at zero theta in the continuum convention with electric energy (g^2/2)integral Phi^2. The bosonized normalization and interval boundary map are given in [Okuda, Secs. II.2 and III.2](https://arxiv.org/html/2210.00297v2). Define a specified spatial source q(z)=sum_j q_j 1_{I_j}(z). Then

    H_C[q]-H_C[0] = g^2 integral q(z)Phi(z) dz + (g^2/2)integral q(z)^2 dz,
    (g^2 a)^(-1) dH_C/dq_j at q=0 = (1/a)integral_{I_j}Phi(z) dz.

The plus-sign source convention is explicitly a background electric-flux shift; translating it to a theta-angle sign must follow the chosen canonical conventions. These derivatives identify the *conjugate observables* for the same prescribed cell source. Equality of continuum and lattice generating functionals or low-energy matrix elements is still a matching claim. It does not follow merely because two derivatives have parallel definitions.

## 3. Realizing the historical B2 instrument

Let F=sum_j w_j y_j+c_B, with c_B=(b_left+b_right)/4 and the unchanged B2 geometry weights. Couple a qubit through P=|1><1|. Choose a pulse u(t) of integral kappa=2 and link flux

    q_j(t) = -u(t) w_j P/(g^2 a).

In the y-source description this yields -u P sum_j w_j y_j plus the known self-energy u^2 P sum_j w_j^2/(2g^2 a). Add the apparatus term

    H_app(t) = -u(t)c_B P - u(t)^2 P sum_j w_j^2/(2g^2 a).

The total interaction is exactly -u(t) P F. If implemented with the conventional L-source, also subtract g^2 a sum_j epsilon_j q_j/4. All these terms are proportional to the detector projector and can change its relative phase; they are not ignorable global phases of the joint system.

In the ideal impulsive limit, the controlled unitary is (I-P) tensor I + P tensor exp(i kappa F). Preparing the qubit in |+> and measuring X produces the historical K_+=(I+exp(i kappa F))/2 and K_-=(I-exp(i kappa F))/2. This supplies an explicit, normalized physical source model for that ideal instrument. A finite-duration pulse also evolves H_a, which does not commute with F, and requires a separate switching-error analysis. D1 does not claim a finite pulse exactly reproduces an impulsive kick.

## 4. Boundaries and finite-spacing corrections

D1 preserves the historical geometry: x_i=(i+1)a, L=(N+1)a, I_j=[x_j,x_{j+1}]. Exact algebra gives y_{-1}=y_{N-1}=1/4 in the neutral sector. With this coordinate assignment the outer flux regions are [0,a] and [Na,L]. Their fixed electric energies add only constants. Algebra alone does not uniquely place continuum Dirichlet endpoints at 0,L. Different continuum-length conventions occur in the literature; changing those endpoints is a different matching prescription, not an exact consequence of the finite chain. D1 fixes the original endpoints and fixed outer-cell values before testing.

Nor does source coupling uniquely choose cell averaging instead of point sampling or another source interpolation. For a smooth scalar profile at a cell center z_j,

    average_{I_j} Phi = Phi(z_j) + a^2 Phi''(z_j)/24 + O(a^4).

This elementary expansion is not a quantum operator theorem. More generally, a source coupling to y_j+c a^2 D_a^2 y_j in the interior is gauge invariant and has the same leading smooth continuum target, while changing finite-spacing matrix elements. Its dimensionless coefficient c is not fixed by zero-source dynamics. Boundary/source counterterms can introduce further differences. D1 selects the unimproved source-conjugate y rule, with no fitted improvement coefficient. The exact normalization and staggered subtraction are now operationally fixed; all finite-spacing matching coefficients and the boundary extrapolation have not thereby been derived.

## 5. Backaction identities and independent tests

The exact commutator [F,T_j]=w_j T_j implies W^dagger T_j W=exp(-i kappa w_j)T_j. For the real, zero-current ground state,

    Delta E_lat = sum_j [cos(kappa w_j)-1] <H_hop,j>.

The direct energy expectation in W|Omega> provides a separate implementation check. In the continuum, W translates the scalar momentum by kappa f_a/sqrt(pi), giving

    Delta E_C = kappa^2/(2 pi) integral f_a(z)^2 dz
              = kappa^2/(2 pi a) sum_j w_j^2.

Here f_a=w_j/a on the interior cells and zero elsewhere. The boundary constant affects the detector phase but not this energy. For either measured branch K_s|Omega>, the unnormalized injected energy is Delta E_W/4, so its conditional value is Delta E_W/(4 p_s). This explicitly includes selective backaction. The cosine contains finite-a corrections beyond its quadratic term; other corrections reside in the state's hopping matrix elements. No kinetic continuum approximation is used to replace the measured lattice result.

Independent prospective probes will use source-native weights v_j=integral_{I_j}g(z)dz, not the quadratic reconstruction. The continuum comparator is explicitly sum_j (v_j/a)integral_{I_j}Phi+c_g I, with c_g=(1/4)integral_outer g. This last term is an identity-valued apparatus calibration: a Dirichlet endpoint value does not make the true continuum field constant throughout a finite-width outer cell. New spatial modes, asymmetric boundary-sensitive probes, first/second moments, insertion matrix elements, vacuum cross-covariances, and energy injection test the dictionary on ground and detector-reached states. The protocol fixes all choices before these observables are evaluated. A finite numerical pass supports only those source observables and states. It cannot establish a unique regulator-independent operator identity or uniform control on every branch of every possible detector history.

For the step smearing f_a, the kicked continuum state has finite mean energy but generally logarithmically divergent energy variance: the jumps give mode coefficients of order 1/n, so sum omega_n f_n^2 need not converge. D1 uses the exact Parseval mean energy, not a truncated H^2 estimate. Continuum energy statements here concern quadratic forms on finite-energy states; they do not assert that the kicked state lies in the operator domain of H.

# Matching and interpretive controls

The action, masses and charge convention are those of PROTOCOL.md. These calculations
are internal controls, not an observation or a new universal law.

## 1. Charge, expansion and the complete energy budget

With X/v=(u+iw)/sqrt(2), q=(u,w,S), chi=a q, and P=chi', the conserved charge is

    Q = integral(chi_w P_u - chi_u P_w) d^3x.

The field-space rotation is a convention; the initial orbit uses P_w=-omega chi_u.
The proper charge density is n_Q/v^3 = sqrt(lambda_phi) Q/(a^3 V_com).
The restored homogeneous orbit has F= kappa^(1/3), omega=sqrt(0.1) F,
so nbar_i=kappa sqrt(0.1). Its k=0 S mode becomes unstable at a=kappa^(1/3).
The finite-k crossing instead uses a^2=F^2+k_lattice^2. Validation compares that
mode against a separate ODE; it does not misidentify growth of a canonical norm
with an internal instability of a finite object.

The spatially discretized Hamiltonian is H_conf = integral[P^2/2+G+V_tau] and

    H_conf' = -a a' integral(chi_S^2-a^2).

Here G is the forward-link gradient energy, which is paired with the seven-point
Laplacian. The periodic flux sums to zero. The KDK charge symmetry is exact up to
roundoff; the work quadrature and integration errors are second order in dt.

For physical density and mean isotropic pressure,

    rho = a^-4 < K + G + V_tau >,
    p   = a^-4 < K - G/3 - V_tau >,
    K   = sum(P-a'/a chi)^2/2.

These include ALL the simulated energy. They do not require a decomposition into
walls, attached lumps and propagating waves. A valid wall/receiver split would need
an additional operational definition and controls for its energy and flux.

## 2. Area weights, energy weights and the one-speed example

For a thin, unoccupied wall element of common rest tension sigma and speed v_n,
dE = sigma gamma dA. Let R=E/(sigma A).

If s_A=<v_n^2>_A is area weighted, convexity gives

    R = <(1-v_n^2)^(-1/2)>_A >= (1-s_A)^(-1/2).

There is NO finite upper bound from s_A alone. For s_A=0.456^2, a mixture of rest
walls and speed-0.99 walls with fast area fraction s_A/0.99^2 gives R=2.29179,
whereas the one-speed value is 1.12362. Letting the fast speed approach one gives
arbitrarily large R at this fixed s_A.

If instead s_E=<v_n^2>_E is genuinely energy weighted,

    R = 1/<sqrt(1-v_n^2)>_E,
    (1-s_E)^(-1/2) <= R <= (1-s_E)^(-1).

The first inequality follows from concavity of sqrt; the second follows from
sqrt(1-x)>=1-x. At s_E=0.456^2 the interval is 1.12362--1.26252. The one-speed
control therefore supplies the LOWER edge, not an upper limit of 12 percent.

For isolated, unoccupied, planar profiles, the ratio of integrals Sdot^2/|grad S|^2
measures this energy-weighted speed squared. It has not been established to do so
for a charged network containing restored cores and radiated scalar waves.
The current suite retains that ratio as a PROXY; values above one diagnose failure
of the thin-wall interpretation, not superluminal propagation.

The smooth coarea kernel integrates to one:

    integral[-1,1] 3(1-S^2)/4 dS = 1.

For a phi^4 planar kink the rest energy measure is proportional to (1-S^2)dS.
The sharp |S|<0.5 energy fraction is 11/16=0.6875; the fraction weighted by
(1-S^2) on |S|<1 is 4/5=0.8. These calibrate the two phi-energy PROXIES on
unoccupied walls. They do not determine the occupied-wall tension or the energy
in X, and cannot be treated as a full wall/receiver partition.

## 3. Pressure and the radiation-equivalent comparison

The exact instantaneous scalar expressions above give

    p-rho/3 = (2K-4V_tau)/(3a^4).

In particular take a homogeneous field exactly at a potential minimum and give it
arbitrarily small nonzero velocity. At that instant rho>0, p=rho, and w=1. Thus
an arbitrary-small MEAN density does not prove an instantaneous w<=1/3 theorem.
Rapid, weakly excited massive oscillations near the minimum have coarse-grained
w approximately zero; relativistic quanta have w approximately 1/3. Those state
and averaging assumptions can support a cosmological argument, but they must be
specified. The counterexample does not demonstrate a long-lived stiff phase or
rescue the persistent-network benchmark.

Similarly Delta N_eff-equivalent is an energy conversion. Using rho_nu,one/rho_gamma
=(7/8)(4/11)^(4/3), rho_hidden=rho_standard_rad corresponds to roughly 7.4 extra
neutrino equivalents. This alone is not a CMB likelihood for walls, matter, a
coupled fluid or an inhomogeneous universe. The stress, perturbations, background
and inference assumptions matter. No such likelihood is computed here.

## 4. Spatial preparation, cross-correlation and the zero mode

For the sparse preparation let xi_u, xi_w, xi_S each have unit variance, with
S=eps_S xi_S, u=F c(1+eps_X xi_u), w=F c eps_X xi_w. The deterministic charge
normalization is c^2=1/(1+2 eps_X^2) for the exact band-limited grid average.
The leading spatially uniform S acceleration from the interaction is

    <chi_S''> = -2 F^2 c^2 eps_S eps_X <xi_S xi_u> + cubic terms.

The spatial mean Laplacian vanishes; zero mean S alone does not remove the
cross-correlation source. For the real Fourier polynomial, the covariance is
one-half the sum of products of its cosine/sine coefficients, which is calculated
in RESULTS.json. This term can seed a coherent zero mode; its sign need not fix
the eventual vacuum because the mode evolves before the instability.

Our 64 mode coefficients are global preparation variables. Keeping their number
fixed as the box grows is not the thermodynamic Gaussian ensemble. The additional
dense-spectrum control addresses this difference. Neither a zero one-point mean,
an ensemble Z2 symmetry, nor a three-seed outcome alone fixes the spatial
correlations or an infinite-volume percolation probability.

## 5. What the extrapolated bias threshold means

Claude's fitted gamma and its extension to growth K over 11--17 decades were known
before this calculation. No run here measures a systematic bias near 10^-23 or
validates that extrapolation. Even within the proposed law eps proportional to
K^(-gamma), log10(eps) changes by -log10(K) delta_gamma. At log10(K)=16.7,
delta_gamma=0.08 alone moves eps by a factor about 22. This does not include
model-form error, a charged-network change of exponent, history or correlations.

Primary comparison read at abstract level: Correia, Leite and Martins,
arXiv:1804.10761 / Phys. Rev. D 97, 083521. It explicitly distinguishes analytic
decay laws for a potential bias from those for biased initial conditions. That
comparison neither proves nor refutes Claude's fitted response exponent.

## 6. Physical matching remains external

At crossing the simulated H/m_X is about 0.01575 (kappa=2, tau_i=40). The quoted
reference histories give about 5.5e-10 to 4.7e-9. This is a dimensionless mismatch;
changing units cannot remove it while keeping v and the yield/history fixed.
The experiments are accelerated classical preparations with a prescribed radiation
background. D=0 is invariant classically, not a test of perturbations in D.

The cold control supplies charge as initial data at the start of the calculation.
It does not construct a later-time source within an action with exact U(1)_X.
The physical scale, net yield, thermal history and vacuum-energy matching remain
external inputs. No new portal, inflaton or charge-producing interaction is adopted.

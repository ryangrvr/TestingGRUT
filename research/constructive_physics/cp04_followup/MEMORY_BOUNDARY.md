# Coupled matching and vacuum-memory boundary

2026-10-11; follows user review of 1b4945a. Analytical follow-up only. No new action, sigma replacement, CP-05, Stage-3 card, or frozen-packet rerun.

## The matching Jacobian is the exact object

Let x=ln lambda, u=ln|v| and A=sigma''/sigma'. On a specified smooth family of local solutions before final global matching, define

F(x,u)=4x+ln(mu^4 Omega(x,u))-ln sigma'(z(x,u)).

Then F=0 selects the global solution and, where Delta=F_x is nonzero,

dx/du=(A z_u-(ln Omega)_u)/Delta,
Delta=4+(ln Omega)_x-A z_x.

Equivalently Delta=partial_x ln[mu^4 Omega_tilde/sigma'(z)], where Omega_tilde=lambda^4 Omega. This derivative is transverse to the final matching surface, on the prior local-solution family. Evaluating it along a fully matched family where the ratio is identically one would instead give zero and lose the selection test.

For a pure constant-v shift at fixed lambda, fixed matched non-vacuum coefficients and v-independent boundary/preparation data, the reduced gravitational/matter equations lose v. On the same locally unique solution branch, Omega_u=0 and historic tau_u=0; z_u=-v/mu^4. Cancellation does not by itself guarantee branch uniqueness or invariant boundary data. Under those qualifications the numerator simplifies, but z_x generally remains.

Write C=historic_tau/(4 lambda^4 mu^4)=Lambda_res/(lambda^4 mu^4). At fixed v,mu,

z_x=C_x=(partial_x Lambda_res-4 Lambda_res)/(lambda^4 mu^4).

**Small C is not a bound on C_x.** For example C(x)=epsilon+L(x-x0) has C(x0)=epsilon and C_x=L; smoothness alone places no bound on L. This is a logical counterexample to an inference, not a cosmological solution. To neglect A z_x physically one must establish

|A C_x| << |4+beta|,

not merely estimate |C|. One needs a response bound for Lambda_res, a bound on A, and control of proximity to beta=-4. Planckian mu alone is insufficient if lambda is small. Also identifying a historic residual with the currently observed dark-energy density requires the model's full matter/cosmology decomposition; it is not automatic. No statement that the correction is large is established either. Computing beta may be sufficient in a demonstrated regular regime; it is not sufficient universally.

## Fluid memory is not yet vacuum memory

The exponential shear kernel eta0 exp(-s/tau)/tau is causal and passive, and the supplied medium equations lead to c_shear^2=eta0/(w tau), w=rho+p>0. This result needs a physical rest frame and positive inertial enthalpy. A Poincare-invariant vacuum expectation T_munu=-rho g_munu instead gives p=-rho and w=0. The same medium equations then have no inertial velocity equation of that form; dividing by w is invalid. This does not imply vacuum correlations have no causal memory. It rules out identifying this specific hydrodynamic speed bound with vacuum gravity without additional physics.

The zero-momentum function eta0/(1-i omega tau) is not a complete Lorentz-covariant spacetime response. It may be a medium response or an effective observer/worldline response. For a vacuum scalar channel, a full retarded correlator depends on invariant momentum with the retarded boundary prescription; tensor stress channels also require conservation identities, projectors and contact terms. Physical spectral weights may contain particle poles and continuum thresholds. Positive frequency dissipation does not require a nonnegative time-domain kernel, so the positive-kernel first-moment tau_eff is a restricted parametrization, not a universal vacuum observable.

A gravitational response also cannot be inferred merely by relabelling a stress correlator: the metric constraints, physical source/readout map and gravitational matching are separate obligations. No bath, preferred velocity, microscopic clock, relaxation scale or generated graviton is introduced here.

## Verified literature boundary

Eling–Guedens–Jacobson, [gr-qc/0602001](https://arxiv.org/abs/gr-qc/0602001), and Chirco–Liberati, [0909.4194](https://arxiv.org/abs/0909.4194), already develop nonequilibrium spacetime thermodynamics/dissipation. These primary abstracts were checked; no claim of independently reproducing their complete derivations.

For vacuum spectral rather than fluid response, [Floerchinger QFT lectures, lecture 25](https://www.tpi.uni-jena.de/~floerchinger/qft2/lecture25/) is an author-hosted teaching source for the vacuum spectral representation; the statement here is ordinary QFT, not a GRUT result. Exact p=-rho follows directly from the supplied Lorentz-invariant stress form, without adopting particular claims about divergence orders in vacuum-energy papers.

## Disposition

Retain the matching Jacobian, reject the unproved universal small-derivative estimate, and separate medium, horizon and physical vacuum response. This refines what a future authorized microscopic construction must calculate. It neither establishes a new theory nor authorizes one. Prior cards and CP-04 remain unchanged; independent review outstanding.

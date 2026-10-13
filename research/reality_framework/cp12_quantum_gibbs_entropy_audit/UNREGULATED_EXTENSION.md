# CP-12: regulator-free Gibbs rarity theorem for narrow spatial graph families

**2026-10-13. ADDITIVE CANDIDATE PROOF, INDEPENDENT PARTNER AUDIT PENDING.**
No changes to frozen CP-12, GRUT-RAI, or Stage-3 Card 3. This is stronger than the finite-cutoff result in REPORT.md: the U(1) electric fields are **unbounded integers**, as in the actual CP-12 law. It addresses only the exact canonical Gibbs state at fixed beta>0 and a graph class satisfying the cardinality bound stated below.

## Conditional theorem

Fix full physical CP-12 in Q=nu N (nu>1), M=N(N-1)/2, with fixed beta>0 and frozen positive mass/electric couplings. For any set S_N of labelled Q-edge graphs with |S_N| <= N! exp(aN) (constant a), let Pi_S retain *all* valid physical fermion/electric states on those graphs. Then

\[
\boxed{\Pr_{\beta}(\mathcal S_N)=\operatorname{Tr}(\rho_{\beta}\Pi_{\mathcal S_N})\to0}
\]

with a quantitative bound

\[
\boxed{\Pr_{\beta}(\mathcal S_N)\le O(1/\log N)+
\exp[-(\nu-1)N\log N+O(N\log\log N)].}
\]

The limit uses *unbounded* U(1) electric fields. The count assumption is satisfied for all relabellings of any one periodic square/cubic lattice (and at most exp(O(N)) unlabelled local variants), **not shown** for a general operationally defined spatial phase.

## Proof: bounded quantum perturbation plus relative entropy

Write H=A+B in the exact physical Hilbert space, with

\[
A=\kappa Q+2m\sum_i n_i+\mathcal E,\quad
\mathcal E=\frac{\epsilon_E}{2}\sum_{e,i\in e} E_{ie}^2,\quad
B=t H_{\rm hop}+H_{\rm rewire}.
\]

Each of the 2Q incidence hopping operators plus adjoint has norm <=1, and the positive rewiring graph Laplacian has maximum weighted degree Gamma Q(M-Q)/M. Therefore

\[
\|B\|\le 2|t|Q+2\Gamma Q(M-Q)/M\le cQ,
\qquad c=2(|t|+\Gamma).
\]

At any finite N, A has trace-class Gibbs exponential since the electric energy grows quadratically in integer fields, and bounded B preserves trace class. Let rho=e^{-beta H}/Z_H and sigma=e^{-beta A}/Z_A. Min-max eigenvalue bounds give |log Z_H-log Z_A|<=beta||B||. Exactly,

\[
D(\rho\|\sigma)=-\beta\operatorname{Tr}(\rho B)+\log Z_A-\log Z_H
\le2\beta cQ.
\]

For K>0 define P_> = 1_{electric energy > KN}; let p=Tr(sigma P_>) and q=Tr(rho P_>). The neutral-sea states on each graph provide Z_A >= binom(M,Q)e^{-beta kappa Q}. Because absent edges carry no rotors, Gauss constraints only eliminate basis states. With Theta(s)=sum_{z in Z} exp(-s z^2),

\[
p\le\binom{N+Q}{Q}\Theta(\beta\epsilon_E/4)^{2Q}e^{-\beta KN/2}
\le e^{-[\beta K/2-c_0]N},
\]

where c_0 is N-independent for Q=nu N. Data processing under the two-outcome P_> measurement gives

\[
D(\rho\|\sigma)\ge d_{\rm binary}(q\|p)
\ge q\log(1/p)-\log2.
\]

Choose K_N=log N; eventually beta K_N/2>c_0. Then

\[
q\le\frac{2\beta cQ+\log2}{(\beta K_N/2-c_0)N}
=O(1/\log N).
\]

Crucially this holds even though P_> does not commute with the full H.

For the complementary low-electric-energy subspace, sum E^2<=2Nlog N/epsilon_E. Gaussian Chernoff counting at eta=1/log N gives

\[
\#\{\mathbf E:\mathcal E\le N\log N\}
\le\exp(2N/\epsilon_E)\Theta(1/\log N)^{2Q}
=\exp[Q\log\log N+O(N)],
\]

using Theta(s)<=1+sqrt(pi/s). Thus

\[
\operatorname{rank}(\Pi_S P_\le)
\le |\mathcal S_N|\binom{N+Q}{Q}
\exp[Q\log\log N+O(N)].
\]

Meanwhile the mutually orthogonal neutral-sea states have exact full-H diagonal energy kappa Q+Gamma Q(M-Q)/M, so spectral Jensen gives

\[
Z_H\ge\binom M Q\exp[-\beta(\kappa Q+\Gamma Q(M-Q)/M)].
\]

Because H_rewire and mass/electric terms are nonnegative and ||tH_hop||<=2|t|Q,

\[
\|\rho\|\le e^{\beta(2|t|-\kappa)Q}/Z_H.
\]

Pi_S commutes with the electric low/high projectors (though none need commute with H). Hence

\[
\operatorname{Tr}(\rho\Pi_S)\le
\operatorname{rank}(\Pi_S P_\le)\|\rho\|+q.
\]

Use log binom(M,Q)=nu Nlog N+O(N), log binom(N+Q,Q)=O(N) and log |S_N|<=Nlog N+O(N) to obtain the theorem. This retains the *complete quantum Gibbs state* without assuming graph-diagonal H.

## Numerical and methodological controls

The conversation research packet includes the standalone stdlib \`unregulated_bound.py\`, which checks theta-function and low/high-energy asymptotic bounds, and \`finite_entropy_validation.py\`, which independently constructs complete small Gauss-law-valid CP-12 matrices and checks the relative-entropy identity/bound and binary data-processing at (N,Q,L)=(3,2,1),(4,2,1),(4,3,1). Markers:

- \`CP12_UNREGULATED_RELATIVE_ENTROPY_BOUND_ALGEBRA_CHECK_PASS\`
- \`CP12_ENTROPY_DATA_PROCESSING_FINITE_H_CHECK_PASS\`

The small finite matrices validate algebraic steps, **not** the unbounded-rotor thermodynamic limit, whose proof relies on trace-class and Gaussian estimates.

## Scope / falsification requirements

1. Only **low-cardinality graph targets** |S_N|<=N! e^{O(N)} are excluded. The result cannot be advertised as an absence of *all* extended or spatially organized graph phases; a broad operationally defined correlated phase could evade it.
2. The proof assumes fixed positive temperature beta>0. Zero-temperature, beta_N growing, nonequilibrium states and cosmological preparation remain open.
3. The CP-12 graph basis is **labelled**. If event permutations are physical gauge redundancies, the quotient Hilbert space and equilibrium ensemble must be rebuilt, not inferred from symmetry alone.
4. CP-12 still conserves Q and has not selected nu. No physical 3D sector, graviton mode, Einstein energy constraint, Lorentz limit or novel empirical prediction is derived.
5. No novelty claim is made for the relative-entropy mathematics or the graph-entropy theme; related prior work: Konopka 2008, PRD 78 044032; Wilkinson & Greentree 2015, PRD 92 084007.

**Partner-audit request:** verify rigorously the uniform boundedness of H_rewire in the physical rotor Hilbert space, trace-class A, the quantum relative-entropy data processing step and physical-Gauss factor counts. Then attempt to define a **broad operationally 3D projector**, count its microstates, and test causal responses in that same fixed-Q Gibbs or justified non-Gibbs state. This is the necessary escape test; do not tune the frozen Hamiltonian.

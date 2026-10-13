# Unregulated CP-12: Gibbs rarity and the resource needed for unitary formation

13 October 2026. Independent analytical follow-up to the [candidate Gibbs audit at 4bb0fa5](https://github.com/ryangrvr/TestingGRUT/blob/4bb0fa5a2835a698168230e5000f31ca0cf880ca/research/reality_framework/cp12_quantum_gibbs_entropy_audit/UNREGULATED_EXTENSION.md). These arguments concern the full labelled physical Hilbert space, unbounded integer rotors, positive mass/electric couplings, and the unchanged Hamiltonian. They are not general no-go theorems for emergent spacetime. Independent mathematical adjudication remains appropriate.

## 1. Reference partition function and bounded interactions

Write

\[
H=A+B,\quad A=\kappa Q+2mN_v+\mathcal E,
\quad\mathcal E=\frac{\epsilon_E}{2}\sum E_{ie}^2,
\quad B=H_{\rm hop}+H_{\rm rewire}.
\]

Let \(G_N=\binom MQ\), \(h_N=2|t|Q\), and \(g_N=\Gamma Q(M-Q)/M\). The norm-one incidence hopping terms give \(\|H_{\rm hop}\|\le h_N\). The physical rewiring operator is a direct sum of finite graph Laplacians at fixed fermion/electric configuration. Every row has weighted degree at most \(g_N\); consequently

\[
H_{\rm rewire}\ge0,\quad \|H_{\rm rewire}\|\le2g_N,
\quad \|B\|\le b_N:=h_N+2g_N=O(N)
\]

in a sparse extensive sector. Unbounded electric fields do not make this rewiring norm unbounded. They remove transitions or change the finite fiber, without increasing its degree. Positivity also gives the sharper form inequality \(A\le H+h_NI\).

For any \(\alpha>0\), define

\[
F_N(\alpha)=\binom{N+Q}{Q}
\left[1+\sqrt{\frac{2\pi}{\alpha\epsilon_E}}\right]^{2Q}.
\]

There are at most \(\binom{N+Q}{Q}\) active fermion choices on each graph. Summing over all \(2Q\) integer electric legs overcounts the Gauss-constrained states, which is safe. The bound \(\sum_{z\in\mathbb Z}e^{-az^2}\le1+\sqrt{\pi/a}\) follows from a monotone Gaussian integral. Each graph has a neutral sea with \(A=\kappa Q\). Therefore

\[
G_Ne^{-\alpha\kappa Q}\le Z_A(\alpha)
\le G_Ne^{-\alpha\kappa Q}F_N(\alpha).
\]

This proves trace class at every finite \(N\). Bounded perturbation preserves a lower-bounded compact-resolvent Hamiltonian and its trace-class Gibbs exponential. Entropy and energy expectations are finite at positive temperature; for example the trace at half the inverse temperature controls the extra energy factor.

## 2. A shorter regulator-free Gibbs proof

Let \(\rho_\beta=e^{-\beta H}/Z_H\) and \(\sigma_\beta=e^{-\beta A}/Z_A\). Bounded perturbation and min-max eigenvalue bounds imply

\[
|\log Z_H-\log Z_A|\le\beta b_N,
\qquad D(\rho_\beta\|\sigma_\beta)\le2\beta b_N.
\]

For a graph set \(\mathcal S\), let \(\Pi_\mathcal S\) include all its physical rotor/fermion states. Its reference probability satisfies

\[
q_\mathcal S=\operatorname{Tr}(\sigma_\beta\Pi_\mathcal S)
\le\frac{|\mathcal S|}{G_N}F_N(\beta).
\]

Set \(L_N=\log(G_N/|\mathcal S|)-\log F_N(\beta)\). If \(L_N>0\), binary measurement data processing and \(d(p\|q)\ge p\log(1/q)-\log2\) give

\[
\boxed{\operatorname{Tr}(\rho_\beta\Pi_\mathcal S)
\le\min\left\{1,\frac{2\beta b_N+\log2}{L_N}\right\}.}
\]

For \(Q=\nu N\), fixed \(\beta>0\), \(\nu>1\), and \(|\mathcal S_N|\le N!e^{aN}\), the denominator is \((\nu-1)N\log N-O(N)\). Thus the probability is \(O(1/\log N)\). This reaches the audit's vanishing-probability conclusion directly, without an intermediate electric cutoff or high-energy split.

It is important not to replace this argument with an operator inequality \(\rho_\beta\le e^{2\beta b_N}\sigma_\beta\). That inequality is false for general noncommuting \(A,B\). `entropy_control.py` includes a two-level countercontrol, as well as actual CP-12 matrices.

## 3. The electric split extends the cooling window

The original audit's split can also be strengthened. The neutral mixture is a Gibbs variational trial state with entropy \(\log G_N\), energy \(\kappa Q+g_N\), and free energy \(\kappa Q+g_N-\beta^{-1}\log G_N\). For the exact Gibbs state,

\[
S(\rho_\beta)\le \frac\beta2\langle A\rangle+\log Z_A(\beta/2),
\qquad \langle H\rangle\ge\langle A\rangle-h_N.
\]

Combining the free-energy inequalities yields

\[
\boxed{\langle\mathcal E\rangle_\beta
\le 2g_N+2h_N+\frac2\beta\log F_N(\beta/2).}
\]

For \(\beta\ge\beta_{\min}>0\), this is \(O(N)\), uniformly also when \(\beta\) grows. Markov's inequality bounds the electric tail above \(KN\) by this expectation divided by \(KN\). The number of states below that threshold on each graph is at most

\[
R_N(K)=\binom{N+Q}{Q}e^{2N/\epsilon_E}
\left[1+\sqrt{\pi K}\right]^{2Q}.
\]

Spectral Jensen on the \(G_N\) orthogonal neutral seas and \(H\ge(\kappa-2|t|)QI\) give

\[
\Pr_\beta(\mathcal S)\le
\frac{2g_N+2h_N+2\beta^{-1}\log F_N(\beta/2)}{KN}
+\frac{|\mathcal S|}{G_N}R_N(K)e^{\beta(h_N+g_N)}.
\]

With \(K=\log N\), \(\beta_N=b\log N+o(\log N)\), and \(\beta_N\ge\beta_{\min}\), the second term has leading exponent

\[
[-(\nu-1)+b\nu(2|t|+\Gamma)]N\log N+O(N\log\log N).
\]

Therefore even the unregulated Gibbs state obeys the sufficient suppression window

\[
b<\frac{\nu-1}{\nu(2|t|+\Gamma)}.
\]

These constants are sufficient bounds, not phase transitions or a claim that cooling faster selects space.

## 4. An all-time nonequilibrium obstruction

Consider the completely unorganized neutral mixture

\[
\rho_{\rm mix}(0)=G_N^{-1}\sum_G|\Omega_G\rangle\langle\Omega_G|.
\]

Under isolated full CP-12 unitary evolution, \(\|\rho_{\rm mix}(T)\|=G_N^{-1}\) and its energy remains \(\kappa Q+g_N\). The form inequality above implies \(\langle\mathcal E\rangle_T\le g_N+h_N\) for **all** times. Since the graph and electric projectors commute,

\[
\boxed{\Pr_T(\mathcal S)\le
\frac{g_N+h_N}{KN}+
\frac{|\mathcal S|}{G_N}R_N(K).}
\]

No temperature, equilibration, or graph-diagonal approximation is assumed here. At \(K=\log N\), this reproduces the audit's \(O(1/\log N)+\exp[-(\nu-1)N\log N+O(N\log\log N)]\) form, now uniformly in time for this nonequilibrium initial state.

A more general entropy condition follows from the same reference state. For any initial state of finite energy \(E_0\) and finite entropy \(S_0\),

\[
D(\rho(T)\|\sigma_\beta)
\le\beta(E_0-\kappa Q+h_N)+\log G_N-S_0+\log F_N(\beta).
\]

If \(E_0-\kappa Q=O(N)\) and \(S_0\ge\log G_N-O(N)\), binary data processing excludes the same narrow targets at all times. Unitary generation of off-diagonal coherence cannot remove the initial entropy.

For fixed nonzero target probability and bounded energy density, the stronger norm/rank argument requires initial min-entropy \(-\log\|\rho(0)\|\le\log|\mathcal S_N|+O(N)\). Thus a narrow \(N!e^{O(N)}\) target in \(\nu>1\) needs a min-entropy deficit of at least \((\nu-1)N\log N-O(N)\) relative to the neutral mixture. This is a preparation requirement, not a rule selecting which low-entropy state nature realizes.

## 5. What the coherent experiment actually supplies

The pure unordered sea \(|\Omega\rangle=G_N^{-1/2}\sum_G|\Omega_G\rangle\) has zero entropy and full-H mean energy \(\kappa Q\). Initially its relative entropy to \(\sigma_\beta\) exceeds that of the neutral mixture by exactly \(\log G_N\), although their energy difference is only \(g_N=O(N)\). It supplies a large entropy deficit without supplying a spatial ordering. The finite formation experiment uses precisely this resource; it has not derived its physical preparation.

The experiment has \(Q=N-1\), whereas the \(N\log N\) exclusion requires \(\nu>1\). These regimes must not be combined as one successful phase. Nor does graph-record dephasing of five already restricted spatial graphs produce the high-entropy mixture over every possible graph.

For a high-entropy, bounded-energy preparation to give a broad graph-defined phase nonvanishing probability, that phase must contain at least \(G_Ne^{-O(N)}\) graphs. For a genuinely quantum operational projector, bare graph cardinality may be inappropriate: the general diagnostic is its probability in \(\sigma_\beta\). Exponentially small reference probability and an \(O(N)\) relative-entropy budget remain incompatible with nonvanishing physical probability.

One exact quantum-code control makes this distinction concrete. Since \(A\ge\kappa Q\) and \(Z_A\ge G_Ne^{-\beta\kappa Q}\), \(\|\sigma_\beta\|\le1/G_N\). Any finite-rank operational projector, including one supported coherently on every graph, has reference probability at most \(\operatorname{rank}(\Pi)/G_N\). Merely spreading a small code subspace over many graph configurations does not remove its entropy rarity. This statement does not bound an arbitrary infinite-rank macroscopic phase; its reference probability must be calculated instead.

None of these results excludes low-entropy pure preparations, entropy exchange with additional apparatus/environment dynamics, or a sufficiently broad operational spatial phase. Permutations-as-gauge require a new physical ensemble; this proof does not substitute labelled counting for that ensemble.

## 6. A specified gauge quotient can be counted separately

There is one conditional extension worth making explicit. Choose the signed permutation-invariant sector constructed in `invariant_evolution.py`: permutations act simultaneously on vertices, edge fermions, and electric legs, with the neutral-background character chosen so neutral graph-orbit sums are admitted. This is an additional physical interpretation, not a gauge principle derived by CP-12.

Let \(K_N\) be the number of unlabelled \(Q\)-edge graph isomorphism classes. There is one orthogonal neutral orbit state for each class, and

\[
K_N\ge G_N/N!.
\]

Every such state has \(A=\kappa Q\) and full-H diagonal energy at most \(\kappa Q+g_N\): internal neutral swaps have nonpositive off-diagonal matrix elements and reduce that diagonal relative to a definite graph. A representative graph's full matter/electric fiber bounds the number and electric partition weight of its invariant states from above. Hence the previous proofs apply with \(G_N\) replaced by \(K_N\), and with the target count replaced by the number of targeted isomorphism classes.

For at most \(e^{O(N)}\) target templates and \(\nu>1\), the rarity denominator again grows at least as \((\nu-1)N\log N-O(N)\). The analogous unregulated Gibbs and high-entropy unitary obstructions therefore survive in **this specified** gauge-sector ensemble. The neutral orbit mixture is a different initial density matrix from the labelled neutral mixture.

A labelled cardinality bound alone need not bound the number of templates if automorphism groups are huge. For connected graphs of uniformly bounded degree, however, fixing a root and its successive neighbor permutations gives an automorphism bound of the form \(N(d!)^N\). In that case the original \(N!e^{O(N)}\) labelled spatial-family bound does imply at most \(e^{O(N)}\) templates. Alternatively that template count can be an explicit premise.

`entropy_control.py` reconstructs the quotient Gibbs ensembles at four and five events and checks their neutral partition trials and relative-entropy identities. This does not audit every possible interpretation of event relabelling.

## Sources and verification

Quantum measurement data processing is established mathematics: [Lindblad, Commun. Math. Phys. 40, 147–151 (1975)](https://doi.org/10.1007/BF01609396). The separable infinite-dimensional form is covered by [Müller-Hermes and Reeb, Ann. Henri Poincaré 18, 1777–1788 (2017), arXiv:1512.06117](https://arxiv.org/abs/1512.06117). The CP-12-specific bounds above are derived here from the frozen law. No general novelty claim is made.

`entropy_control.py` checks the relative-entropy identity, binary data processing, unitary entropy/energy preservation, and the norm/rank electric split on complete finite CP-12 matrices of dimensions 30, 147, and 688. Those calculations validate algebraic steps; the unbounded thermodynamic statements follow from the estimates above and still warrant independent proof review.

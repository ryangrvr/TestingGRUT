# CP-12: unfiltered relational response and circulation capacity

Parent: `7d414961a56a49df1265b9501cd305331354ffc6`. The frozen Hamiltonian, couplings and supplied neutral preparation are unchanged. These are fixed-finite-system results and explicitly bounded finite-time calculations. None is a large-time Taylor extrapolation.

## 1. Physical observables and the absence of an initial loop filter

There are N vertex modes a_i, M=N(N-1)/2 potential edge modes b_e, Q present graph edges, and integer half-leg electric fields. The constraints are

\[
\sum_{e\ni i}E_{ie}=n_i,\qquad E_{ie}+E_{je}=1-n_{b_e},
\qquad \sum_i n_i+\sum_e n_{b_e}=M.
\]

An absent edge has occupied b_e and zero fields. Consequently the number of holes on present edge modes is n_v=\sum_i n_i. The fixed couplings are kappa=1/2, m=t=1, epsilon_E=1/4 and Gamma=1/10. The Hamiltonian consists of graph cost kappa Q, mass 2m n_v, electric energy (epsilon_E/2)sum E^2, the frozen oriented gauge-covariant hopping, and the neutral all-to-all rewiring Laplacian R.

Let rho_w be the **unfiltered** state after the declared source-free wait. Choose source S=n_0 and detector

\[
D=\frac1{N-1}\sum_{j\ne0}n_j.
\]

Let C=C_{not\,0} be the configuration projector onto the existence of at least one directed electric-circulation cycle avoiding vertex zero. An edge on such a cycle has two opposite-sign, nonzero half-leg fields, consistently directed around the cycle. It may contain an edge-fermion hole, for example fields (2,-1). C depends only on non-source-edge variables. For N=4, C is exactly the circulation event on the fixed triangle (1,2,3). For N=5, it is the union event over cycles on the other four vertices, not a tracked named loop.

The physical response path is rho_w, then exp(-i eta S), then the unchanged full H for lag tau. **There is no initial C or occupied-flux filter.** Define

\[
f(\eta)=\langle CD\rangle_\eta,\qquad g(\eta)=\langle C\rangle_\eta.
\]

The final-readout conditional density response is

\[
\boxed{\chi_{D\mid C}=\left.\frac{d}{d\eta}\frac{f(\eta)}{g(\eta)}\right|_0
=\frac{f'(0)g(0)-f(0)g'(0)}{g(0)^2}.}
\]

This measures response within the circulation readout rather than only a change in its probability. Final conditioning is still a readout postselection. It is neither an initial preparation nor an autonomous selector or observer.

For a projector source S, every expectation has the form a+b cos(eta)+c sin(eta). Thus f'(0)=[f(pi/2)-f(-pi/2)]/2, and likewise g' and the plain density response. This is an exact parameter identity. The derivative of f/g itself must use the quotient rule; taking a phase difference of the two ratios would be incorrect.

## 2. Exact short-lag remote-response operator

For any appropriate finite-moment input state,

\[
\chi_O(\tau)=i\operatorname{Tr}\rho_w[S,e^{iH\tau}Oe^{-iH\tau}].
\]

For the off-source average D the constant and linear lag terms vanish as operator identities. Put

\[
W_{0j}=n^{graph}_{\{0,j\}}a_0^\dagger U_{0e}U_{je}^\dagger a_j,
\qquad
K_0=\frac{i t^2}{2(N-1)}\sum_{j\ne0}(W_{0j}-W_{0j}^\dagger).
\]

Each W is gauge invariant: its two field shifts compensate the changes in the vertex occupations and leave the edge Gauss equation unchanged. The orientation signs in the frozen hopping make the two-step single-particle matrix element -t^2 U_0e U_je^dagger. Direct expansion gives

\[
\boxed{-\frac i2[S,[H,[H,D]]]=K_0,\qquad
\chi_D(\tau)=\tau^2\langle K_0\rangle+O(\tau^3).}
\]

Equivalently the coefficient is -t^2/(N-1) times the sum of the imaginary gauge-dressed intervertex coherences. The mass, electric and rewiring terms make no direct contribution at this order. They still affect the waited state and later dynamics.

The circulation projector commutes with every hopping half-leg on an edge incident to zero, as well as with the source. Its first two source-response coefficients vanish. The first commutator changing C uses modes entirely on edges avoiding zero; a source-zero half-leg cannot meet it in only one further commutator. Therefore

\[
\boxed{\chi_C(\tau)=O(\tau^3),\qquad
-\frac i2[S,[H,[H,CD]]]=C K_0.}
\]

Write alpha=Tr(rho_w C K_0) and p_C=Tr(rho_w C)>0. Then

\[
\boxed{\chi_{D\mid C}(\tau)=\frac{\alpha}{p_C}\tau^2+O(\tau^3).}
\]

This coefficient is an expectation of a bounded Hermitian operator in the original unforced state. No physical initial measurement is required to define or compute it. The main neutral sea has finite Hamiltonian moments; its finite-time evolution supports the derivatives used here. These fixed-system expansions do not control long times or a growing-N limit.

## 3. Exchangeability hides a relational phase correlation

The neutral initial ray and H preserve vertex permutations. Swapping zero and j maps W_0j to its adjoint. Hence every unweighted pair coherence has real expectation, and

\[
\boxed{\langle K_0\rangle=0,\qquad \chi_D(\tau)=O(\tau^3)}
\]

in these permutation-invariant unforced states. C distinguishes a relational role: a loop avoiding the source. It does not preserve the source-detector swap, and alpha need not vanish.

There is no preferred label hidden here. The averaged operator

\[
\overline A=\frac1N\sum_i C_{not\,i}K_i
\]

is permutation invariant and gauge invariant. Its expectation in an exchangeable state equals alpha for any chosen source label. Thus the observed positive alpha is an unforced relational correlation, even though every labelled edge has the same unconditional marginal.

In the graph/occupation/electric configuration basis H is real, while C K_0 is purely imaginary Hermitian. Every real density operator therefore has alpha=0. In particular both existing physical energy references are real and have exactly zero quadratic coefficient. Their finite-lag responses can still be nonzero.

Complex conjugation of any state reverses alpha while preserving all configuration probabilities and the entire energy-measurement distribution: real H has real spectral projections. For the initially real neutral sea the conjugated evolved state is its negative-time counterpart. This is a mathematical conjugation statement, without an additional claim about a derived relativistic time-reversal operation. A positive alpha demonstrates phase–circulation correlation, not a preferred labelled graph, dimension or arrow-of-time law.

## 4. Global coefficient norm and its size limitation

At fixed graph and rotor angles, K_0 is a fermionic star linking vertex zero to d_0 present neighboring vertex modes. The star has only one positive and one negative one-particle eigenvalue, each of magnitude |t|^2 sqrt(d_0)/[2(N-1)]. Its many-fermion norm has the same maximum magnitude. Therefore on the unregulated physical space

\[
\boxed{K_0^2\le\frac{t^4 d_0}{4(N-1)^2},\qquad
\|K_0\|,\|CK_0\|\le\frac{t^2}{2\sqrt{N-1}}.}
\]

C commutes with K_0 and d_0. This proves a global bound, rather than using only a retained-space eigenvalue. If the waited-state norm error is e, the alpha expectation error is at most t^2(2e+e^2)/[2sqrt(N-1)]. These bounds give the saved full-law positive coefficient intervals at N4 and N5.

There is a further all-size constraint. Exchangeability and fixed Q give <d_0>=2Q/N. Cauchy-Schwarz yields

\[
|\alpha|^2\le\frac{t^4 p_C\langle C d_0\rangle}{4(N-1)^2}
\le\frac{t^4 p_C(2Q/N)}{4(N-1)^2}.
\]

At Q=nu N with fixed nu, alpha is O(1/N). If p_C has a positive size-independent lower bound, the conditional quadratic coefficient is also O(1/N). Even for rare readouts, the global operator bound already forces this conditional coefficient to be O(N^-1/2). This average-over-all-other-vertices detector is diluted as the system grows. It cannot be extrapolated to a finite macroscopic signal amplitude. The bound does not exclude finite local-detector responses or an emergent causal region; those require relationally resolved detectors and their own controlled size dependence.

## 5. Recording a candidate loop is a physical change

For I=C or the occupied-flux-cycle event avoiding zero, let a=I psi_w and b=(1-I)psi_w. For each phase path and final positive detector O,

\[
\langle O\rangle_{\rm coherent}
=\langle O\rangle_a+\langle O\rangle_b
+2\operatorname{Re}\langle U_\eta a,O U_\eta b\rangle.
\]

Replacing the coherent state by I rho I+(1-I)rho(1-I) discards the cross term. The same decomposition holds for the phase derivative. The packet propagates both subnormalized branches, preserving their probabilities, rather than comparing arbitrarily renormalized selected states.

At N4, wait10 and lag0.5 the density-response interference is resolved for both record choices. At N5 the corresponding interference intervals include zero. The result is that the proposed record branches are **not exactly a classical mixture** in the N4 diagnostic. It does not exclude an approximate coarse-grained classical description at a declared tolerance or in a different scaling limit.

For these non-source-edge records, [I,K_0]=0. Recording leaves the leading quadratic remote-response coefficient unchanged; its effect on the plain susceptibility begins at O(tau^3). The observed finite-lag interference is consistent with that identity. The within-record phase correlation alpha must not be confused with interference between the two record sectors.

For a pure input the trace distance to the pinched state is sqrt[p(1-p)], and the measured-system pinching entropy is h_2(p). The nonselective mean recording work is

\[
\Delta E_{record}=-2\operatorname{Re}\langle I\psi_w,H(1-I)\psi_w\rangle.
\]

Only hopping contributes because I commutes with R and the diagonal Hamiltonian. Source-kick work also comes only from hopping. A global hopping norm bound |t|sqrt(2NQ) gives the stated conservative work error 2|t|sqrt(2NQ)(2e+e^2). Including an apparatus in a unitary description does not create total entropy merely by recording. Its construction and a locally accessible internal observer remain unprovided.

## 6. Exact fermion-response accounting

Let h_e=1-n_be, including zero on absent edges. Define H_inc as mean h_e on the N-1 edges incident to the source and H_other as mean h_e on the remaining (N-1)(N-2)/2 edges. Pointwise number conservation gives

\[
S+(N-1)D-(N-1)H_{inc}
-\binom{N-1}{2}H_{other}=0.
\]

The same identity holds for source derivatives, even after multiplying all terms by a common final C. The saved plain-response Ward residual is approximately 10^-16. Most of the source's local response is balanced by changes in edge occupancy. A vertex-density susceptibility is a matter response, but does not alone measure migration of a previously tagged vertex particle or a causal speed. Setting hopping to zero after the same preparation makes all these density-phase responses exactly zero. The initial vertex-vacuum sea is likewise inert to this phase source before waiting.

## 7. Electric capacity for a spanning, cycle-rich network

Let mathcal C be the union of all directed electric-circulation cycles in a physical configuration. Denote its numbers of edges, nonisolated vertices and connected components by e_C,v_C,c_C, its hole count by h_C, and its cycle rank by beta_C=e_C-v_C+c_C.

An occupied edge in mathcal C has electric square at least 2. A hole edge in mathcal C has square at least 5, since its two nonzero opposite-sign integer fields sum to 1. Every hole outside mathcal C costs at least one. There are n_v holes in total. Summing disjoint contributions proves the unregulated diagonal inequality

\[
\boxed{\sum E^2\ge2e_C+n_v+2h_C.}
\]

For a connected spanning circulation union,

\[
\boxed{\beta_C\le\frac12\sum E^2-N+1-\frac12n_v-h_C.}
\]

An electric budget B therefore must satisfy B>=2(N-1+beta_C)+n_v+2h_C to even represent that configuration. These are necessary bounds, not sufficient formation criteria or sharp minima for every topology.

Several consequences matter for the next experiment:

* B<2N excludes every spanning circulation configuration from the numerical projection. Its projected zero is not full-law absence.
* At B=2N a spanning circulation sector can only have a neutral unit-flux ring as its circulation support, with n_v=0. The full graph may also contain neutral connections. A density-phase source is then inert in that selected sector. Simply reaching the first spanning budget can therefore introduce a second diagnostic obstruction.
* If every one of Q edges circulates, all n_v holes belong to mathcal C, so sum E^2>=2Q+3n_v. A fully circulating degree-four network has Q=2N and needs B>=4N+3n_v. Degree four is an example of a network capable of three-dimensional growth, not a sufficient or necessary dimension test.
* At Q=N a connected spanning circulation network remains a simple ring. Enlarging that sequence does not test a three-dimensional circulating phase.

Electric energy is (epsilon_E/2)sum E^2. These are component-energy costs; negative hopping can compensate them in total energy, so the inequalities do not forbid autonomous formation. They constrain resources and numerical capacity. The independent checker tests 6,075 raw physical states and 414,568 retained orbit configurations across the four existing cases.

## 8. Finite-time projection errors and independent verification

The numerical projection is the unchanged complete reachable electric-square budget with every outgoing full-law transition in B=(1-P)HP. The signed symmetry quotient is exact. The source-fixing basis adds no Hamiltonian term. Its inclusion intertwines the invariant and source-fixed Hamiltonians.

For each phase path, e_path<=e_wait+sum_segment epsilon_segment, with the inherited fifth-order midpoint upper sum of the Duhamel leakage integral. Matrix exponentials evolve the physical projected state; the midpoint Taylor polynomial only bounds leakage. Mixed-state purifications use Frobenius errors.

Every positive contraction O commuting with P obeys

\[
\langle O\rangle_{full}\in
[(\sqrt{p_{projected}}-e)_+^2,\min(1,(\sqrt{p_{projected}}+e)^2)].
\]

Interval arithmetic on these analytic inequalities bounds f,g,f',g' and their quotient derivative. For two subnormalized record paths with errors e_a,e_b the complex cross term has error at most e_a||O b||+e_b||O a||+e_a e_b. These are arithmetic on ordinary floating-point numbers; there is **no formal interval-roundoff certificate**.

The independent checks include a 6,075-state raw Gauss Hamiltonian, 1,875 interior columns for unregulated short-lag identities, signed orbit isometry/intertwiner/outgoing-Gram checks, 18 raw-versus-quotient response comparisons, and 24 separate dense unitary dilations with pure and mixed inputs. The verifier freshly reconstructs N4 and N5 main-state calculations and distinguishes those reruns from hash-validated larger capacity and reference evolution files.

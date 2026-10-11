# Independent CP-07 / CP-08 review: locality must be a physical response, not a link cap

2026-10-11. Standalone partner review of frozen CP-07 at \`48f9c156bd9b93e24d3fdb11def7ffb37b547225\` and CP-08 at \`e3cd89a16478ed3e36debee4373c7a5f33438108\`. This review branch starts from CP-08 and adds new files in **research/reality_framework/** only. This is not CP-09, a modified generating law, a spin-2 completion, or Stage-3 Card 3; no frozen construction or PR is altered.

## What survived direct inspection

CP-07 defined a normalized matter-conditioned relation-growth circuit. Its phases do not alter final diagonal history probabilities because its relation registers are permanent distinguishable which-history labels. The conditional edge probability has a uniform lower bound 1/4; greedy-path/martingale arguments make its height extensive in N, incompatible with a fixed-dimensional Minkowski Alexandrov interval's N^(1/d) chain scale.

CP-08 replaced permanent path records with a Hermitian Hamiltonian on all labelled finite partial orders tensored with matter. It earns recombination of distinct routes, a nontrivial small-loop unitary holonomy, observable record effects and shared finite-system response/noise from the SAME supplied dynamics. Its exact variational bound E_0(N) <= -J floor(N^2/4)+N|Omega| reveals that a fixed positive hopping scale J cannot furnish a stable extensive thermodynamic energy limit on the declared ensemble. This bound does not prove that physical gravitational energy must be extensive or that the exact global ground state is itself a layered order. Neither construction derives a physical clock, conservation of spacetime stress-energy, a gauge-invariant spin-2 mode, a Lorentzian metric or observed gravitational coupling.

## New diagnostic: do not equate Lorentzian locality with a bounded number of links

Let a four-dimensional Minkowski spacetime be Poisson sprinkled at constant four-volume density rho>0, and condition (Palm convention) on an element at x=0. A future element y is a **link** if the Alexandrov interval from x to y contains no other sprinkled element. For proper-time separation tau,

\[
V_4(x,y)=\frac{\pi}{24}\tau^4,\qquad
P(x\prec^\ast y\mid y)=\exp\left[-\rho\frac{\pi}{24}\tau^4\right].
\]

Introduce future-hyperbolic coordinates

\[
t=\tau\cosh\chi,\quad |\mathbf x|=\tau\sinh\chi,\quad
d^4y=\tau^3\sinh^2\chi\,d\tau\,d\chi\,d\Omega_2.
\]

The expected number of future links, integrated over **all** timelike tau but rapidity limited only as a regulator to 0<=chi<=Xi, is

\[
\begin{aligned}
E_\rho[L(0\le\chi\le\Xi)]
&=\rho 4\pi\int_0^\infty d\tau\,\tau^3
e^{-\rho\pi\tau^4/24}
\int_0^\Xi d\chi\,\sinh^2\chi\\
&=\boxed{6\sinh(2\Xi)-12\Xi}.
\end{aligned}
\]

The tau integral is 6/pi, cancelling rho in the final expression. The formula is exact for this unbounded Minkowski Poisson model, independent of rho. It grows as 3 exp(2Xi) when Xi becomes large, so the unbounded future-link expectation diverges. Xi itself is a noncovariant rapidity regulator, **not** a physical invariant cutoff. This is the familiar Lorentzian noncompact-hyperboloid problem, not a novel GRUT discovery. It does not mean that the finite-dimensional CP-08 Hamiltonian already obeys any Lorentz symmetry. It also must not be generalized to arbitrary background, finite sampling region, or a claimed uniformly drawn link without specifying a regulator.

For comparison, in 1+1 Minkowski where V_2=tau²/2 and the rapidity coordinate ranges -Xi..Xi, the same exact expectation is 2Xi.

The existing Bombelli–Henson–Sorkin result proves, under its Poisson-sprinkling and equivariance assumptions, that one cannot construct an intrinsic finite-valency graph while preserving Lorentz invariance. Therefore **a hard finite degree cap on causal links is not a legitimate general solution to the CP-08 locality obstruction** if the intended continuum correspondence uses Lorentz-invariant sprinklings.

This *does not* imply that observable propagation must be nonlocal at long distances: an operator built from nonlocally distributed links can converge to a local wave operator through cancellations, smearing, or intrinsic distance constructions. Benincasa and Dowker (2010) constructed a continuum-corresponding retarded nonlocal d'Alembertian. Boguñá and Krioukov (PRD, 22 Jan 2026) report a local causal-set d'Alembertian defined via intrinsic timelike and spacelike distance reconstruction, with a continuum limit for arbitrary Minkowski test fields. Those results are **external prior art**, not mechanisms selecting an ensemble. A new October 2026 arXiv preprint on link sampling biases highlights why fixed-window numerical link statistics can be misleading; its abstract has been checked, not its derivations.

## CP-08 transport curvature has a precise support

The CP-08 matter generator on a directed relation is P_ij=Z_i X_j and U_ij=exp(-i gamma P_ij). For two independent relations (i,j) and (k,l), the generators commute unless one relation's target becomes the other's source (j=k or l=i), i.e. unless a composable causal two-step path exists; self-reversing pairs are excluded in a poset. Consequently every bipartite height-two sector with only lower-to-upper arrows has mutually commuting transport operators even if it has O(N²) edges. Its elementary Wilson loop is trivial.

In the composable three-event control P_01=Z0 X1 and P_12=Z1 X2 anticommute. Let A=U_01 and B=U_12. Then

\[
\frac{1}{8}{\rm Re}\,{\rm tr}(A^\dagger B^\dagger A B)=1-2\sin^4\gamma.
\]

At gamma=pi/4 this is 1/2. For a noncomposable height-two pair P_01 and P_02 the same expression equals 1 exactly. These follow from standard Pauli algebra and were checked using independently assembled 8x8 matrices. Therefore **non-Abelian matter curvature in this law detects composable chains but supplies no selection energy over the flat layered hypercubes**. This is an exact clarification of an obstruction already present in the CP-08 report, not another prospective generating law.

## Reproducible independent control

[check_link_intensity_and_holonomy.py](check_link_intensity_and_holonomy.py) evaluates independent radial and rapidity integrals with SciPy, compares against the exact 4D formula for rho=0.1, 1, 10 and Xi=0.5, 1, 2, 3, and checks both commuting and anticommuting Pauli-transport loops. It was executed in the analysis environment before publication and passed. Examples:

| Xi | E[number of future links with rapidity <=Xi] |
|---:|---:|
| 0.5 | 1.0512071619 |
| 1 | 9.7611624471 |
| 2 | 139.7395031828 |
| 3 | 1174.2789442217 |

These are *expectations in an externally supplied Poisson/Minkowski comparison model*, not predictions of CP-08's finite Hamiltonian and not an empirical test of spacetime. Their growth warns against treating finite-volume link counts or bounded degree as a Lorentz-invariant criterion.

Run: \`python check_link_intensity_and_holonomy.py\` (NumPy, SciPy).

## The constructive next test (no new candidate authorized)

The next decisive test should operationally define locality and response:

1. Supply **independently calibrated matter clocks** and two detectors with fixed relational definitions; do not choose their support after inspecting an emergent geometry.
2. For a prepared state and one finite evolution law, derive the retarded cross-response \(\chi_{AB}\propto i\theta(t)\langle[A(t),B(0)]\rangle\), its noise, and record correlations. Evaluate the same operators on independently identified continuumlike and layered families, with consistent physical measures.
3. Require a controlled regime with approximate microcausality relative to an **independently reconstructed** spacetime causal relation and a stable continuum wave-operator/dispersion limit. Avoid a hard causal-link degree cap; density, normalization, boundary/window bias, and boost behavior must be disclosed.
4. Only if propagation is earned should one test whether the joint Hamiltonian and state dynamically favor the continuumlike class over the full nongeometric ensemble. A two-history comparison and a scalar field on a supplied geometry are not enough.
5. Finally recover gravitational first-class constraints, universal spin-2 stress-energy coupling and empirical gravitational predictions; standard KMS and finite-system response identities alone are not novelty.

The 2026 local d'Alembertian is a concrete **existing correspondence benchmark**, not evidence that the CP-08 Hamiltonian already produces it. Do not insert that operator as a new geometry-selection weight without specifying and testing the microscopic law.

**Scientific status:** CP-08 is a real finite quantum model with nontrivial recombination. Its next obstacle is relational, Lorentz-compatible *propagation locality* and class-level selection, not merely the number of graph edges. This review adds one exact continuum-comparator calculation and an independently verified holonomy control, not a fundamental theory.

## Sources checked

- CP-07 [frozen report](https://github.com/ryangrvr/TestingGRUT/blob/48f9c156bd9b93e24d3fdb11def7ffb37b547225/research/constructive_physics/attempt07/REPORT.md).
- CP-08 [frozen report](https://github.com/ryangrvr/TestingGRUT/blob/e3cd89a16478ed3e36debee4373c7a5f33438108/research/constructive_physics/attempt08/REPORT.md).
- Bombelli, Henson, Sorkin, *Discreteness without symmetry breaking: a theorem*, https://arxiv.org/abs/gr-qc/0605006 (2006/2009), statement of finite-valency obstruction checked.
- Benincasa, Dowker, *Scalar Curvature of a Causal Set*, https://doi.org/10.1103/PhysRevLett.104.181301 (2010), published abstract checked.
- Boguñá, Krioukov, *Local d'Alembertian for causal sets*, https://doi.org/10.1103/np89-qzjp (22 Jan 2026), publisher abstract checked; full derivation not replicated.
- Boguñá, Krioukov, *Sampling biases in sprinkled causal sets*, https://arxiv.org/abs/2610.03247 (2 Oct 2026), arXiv bibliographic abstract checked; detailed derivation not replicated.

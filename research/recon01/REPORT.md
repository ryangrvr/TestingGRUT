# Reconnaissance 01: a working discriminator, no new theory established

10 October 2026. Development branch `codex/autonomous-reality-research`.
Parent: `eb32a510c0b5f363bcc53eb4bc737feb950ba7f7`.

**Retained result:** a finite, independently measurable comparison of coherent
interaction and noise excludes a defined non-entangling model class. The proof
below is exact in its stated model. **No new GRUT generating law is established.**
One card has been consumed; two remain. No Card-1 repair or Card-2 evaluation
is concealed in this reconnaissance. See [SCOPE.md](SCOPE.md).

## 1. What the primary-source search actually retained

| Published construction | What works | Supplied structure and unresolved scope | Disposition in this pass |
|---|---|---|---|
| [Bisio et al., free QCA fields](https://arxiv.org/pdf/1601.04832), §§II, V–VII | Conditional Weyl, Dirac and Maxwell limits; strong classification within the declared automaton class | Fermionic cells and their algebra, discrete updates, homogeneity/isotropy, Abelian group and minimal dimension are premises. The Dirac mass parameter remains free. Photon commutation is approximate in a restricted occupation regime. The inspected result does not produce the GRUT identity/readout chain. | Retain as a generative construction conditional on inventory; no adoption as a new GRUT law |
| [Bisio et al., interacting Thirring QCA](https://arxiv.org/pdf/1711.03920), §§I–III | Exact two-particle scattering and bound-state structure | On-site interaction phase and free automaton mass are supplied. This is an interacting extension, so the objection is not that every QCA is free. The inspected papers do not establish a gravity–record restriction with an admissible conventional violating point. | Do not spend a card merely copying its update rule |
| [Oppenheim et al., CQ trade-off](https://arxiv.org/pdf/2203.01982), §III; [Layton et al., weak field](https://arxiv.org/pdf/2307.02557), §§2, 7–8 | Consistent hybrid dynamics and conditional decoherence/diffusion bounds | Classical/quantum split, gravitational coupling, field kernels and Markovian approximation supplied. The weak-field paper explicitly distinguishes non-Markovian effective regimes. Full fundamental general-relativistic consistency is not certified by this audit. | Retain as a quantitative alternative to test, not a discovered selector |
| [Fabiano et al., noise discriminator](https://arxiv.org/html/2603.26075v2), §II.2 and B.2 | Coupling/noise comparison including correlated noise | Positive-rate GKSL dynamics and the measurement reduction are load-bearing. This is established-physics selection under additional model premises. | Work through the finite comparator and its measurement certificate |

This is a bounded comparison of these papers, not a survey proving every
possible theory fails. A residual coupling is disqualifying only for observables
it changes; it need not prevent other predictions.

## 2. Exact finite-model statement

**POSTULATED standard comparison class:** two labelled qubits,
`Z_A=Z tensor I`, `Z_B=I tensor Z`; an external clock; Born probabilities;
supplied preparation and calibrated Pauli readouts. All coefficients below have
units of inverse seconds. No structure is called generated.

\[
\dot\rho=-ig[Z_AZ_B,\rho]+\sum_{i,j=A,B}C_{ij}
\left(Z_i\rho Z_j-\tfrac12\{Z_jZ_i,\rho\}\right),
\quad C=\begin{pmatrix}a&c\\c&b\end{pmatrix},
\quad a,b\ge0,\ ab\ge c^2.
\]

The last condition makes the ordinary evolution completely positive. It does
**not** imply that the channel cannot generate entanglement. The generators
commute with the supplied diagonal Hamiltonian; populations, trace and its
energy expectation are conserved. This finite reduction is not a complete
relativistic or thermal gravity model; KMS/FDT for an unspecified bath is not
claimed. Its noise can have an ordinary environment implementation.

**THEOREM (standard-model discriminator).** Within this class, evolution
preserves separability of every initially separable two-qubit state for every
\(t\ge0\) if and only if

\[
\boxed{ab-c^2\ge g^2.}
\]

**Proof of sufficiency.** Let \(\Theta_B\) be transposition of the second qubit
in its Z basis. Direct substitution on every matrix unit gives

\[
\Theta_B\mathcal L\Theta_B=\mathcal D_{C^{\rm PT}},
\qquad C^{\rm PT}=\begin{pmatrix}a&-c-ig\\-c+ig&b\end{pmatrix}.
\]

Under the boxed condition this is a positive-rate GKSL generator, hence
\(e^{t\mathcal D_{C^{\rm PT}}}\) is CP. A separable input has a positive
partial transpose, which stays positive under this evolution. For two qubits,
positive partial transpose is equivalent to separability.

**Proof of necessity.** Start from \(|++\rangle\). Compress the derivative of
the partially transposed state to its zero-eigenvalue subspace, in the ordered
basis \(|-+\rangle,|+-\rangle,|--\rangle\). The result is
\(C^{\rm PT}\oplus0\). Its smaller eigenvalue is

\[
\mu_-={a+b-\sqrt{(a-b)^2+4(c^2+g^2)}\over2}.
\]

If the boxed condition fails, \(\mu_-<0\), giving negative partial transpose
at sufficiently small positive time. In particular the exact right derivative
of negativity is

\[
\dot{\mathcal N}(0^+)=\max\left(0,
{\sqrt{(a-b)^2+4(c^2+g^2)}-a-b\over2}\right).
\]

This proves an entire family statement, not an inference from sampled time
points. Partial-transpose positivity is the standard two-qubit criterion
([Peres](https://arxiv.org/abs/quant-ph/9604005),
[Horodecki et al.](https://arxiv.org/abs/quant-ph/9605038)). Force/noise witnesses
predate this audit ([Kafri–Taylor](https://arxiv.org/pdf/1311.4558),
[Kafri–Taylor–Milburn](https://arxiv.org/pdf/1401.0946)). The finite-class calculation
agrees with the 2026 paper's B.2 after coefficient matching. This is no GRUT
novelty claim.

**Hostile conventional point:** \(g=1,a=b=3/4,c=0\) is CP, energy conserving
in this finite model, and violates the non-entangling inequality. Its initial
negativity rate is exactly \(1/4\,\mathrm{s}^{-1}\). Thus the restriction is
nondefinitional and excludes a physically admissible ordinary quantum channel.
The premise buying it is *non-entangling dynamics*, not a new universal law.
At \(a=b=5/4,c=3/4,g=1\), the condition is saturated. At \(a=b=1,c=3/4,g=1\),
ordinary CP still holds but the interaction entangles. Correlations cannot be
discarded from the precise comparison.

## 3. Independent operational observables and a finite decision

The following repertoire uses the same clock, two devices and calibrated
tomography for every comparison:

1. Ramsey preparation on A with B in each Z eigenstate. The conditional phase
   frequency difference has magnitude \(4|g|\); the normalized contrast envelope
   decays at \(r_A=2a\). Repeat with the roles reversed to measure \(r_B=2b\).
2. Prepare the product state \(|++\rangle\), use joint Pauli tomography, and
   extract the \(00\!\leftrightarrow\!11\) and \(01\!\leftrightarrow\!10\)
   matrix coherences. Their decay rates are
   \(R_+=2(a+b+2c)\) and \(R_-=2(a+b-2c)\), so
   \(c=(R_+-R_-)/8\). These two coherences have no ZZ phase difference.
3. Compare negativity from the same joint state to the rate prediction. Hold
   preparation, calibration and geometry fixed. Collect separated time windows
   and intervention histories to test the assumed semigroup description.

This is a protocol specification, not performed laboratory work. Records,
readout equivalence, memory and subsystem identity are supplied. No common
carrier or \(\epsilon_R\) verdict is earned. Gravity attribution additionally
requires independent bounds on electromagnetic, mechanical and matter-exchange
cross-talk and independently calibrated mass/geometry; a generic ZZ experiment
alone tests a channel rather than gravity.

Given certified boxes \([a_-,a_+]\), \([b_-,b_+]\), \([c_-,c_+]\), \([g_-,g_+]\),
define \(\underline{|c|}\) and \(\underline{|g|}\) as their smallest absolute
values. With a separately certified squared-rate model error \(E\ge0\),

\[
\underline{|g|}^2+\underline{|c|}^2-a_+b_+-E>0
\]

excludes the specified non-entangling class conservatively. Inconsistent CP
calibration is rejected first. Failure to exclude is inconclusive. The code
demonstrates these decisions on exact rational **synthetic** boxes; it does not
manufacture physical error certification from a fit. Model applicability,
finite-time derivative error and noise attribution remain separate obligations.

## 4. Numerical benchmark audit and geometry error

For two equal point masses, mean separation \(d\), superposition displacement
\(s<d\), the exact ZZ coefficient of the four Newtonian branch energies is

\[
g_{\rm point}={Gm^2s^2\over2\hbar d(d^2-s^2)},
\qquad g_{\rm quadratic}={Gm^2s^2\over2\hbar d^3}.
\]

This follows from \((E_{++}-E_{+-}-E_{-+}+E_{--})/(4\hbar)\).
The relative point-geometry truncation error is exactly
\((s/d)^2/[1-(s/d)^2]\), or \(1/99\) at \(s/d=0.1\).
This controls one approximation only, not wave packets, finite radii, shielding
or environmental dynamics.

Using the explicitly rounded [CODATA SI literals](https://physics.nist.gov/cuu/pdf/all.pdf)
in audit.py, with \(s=100\,\mu m\)
and \(d=1\,mm\), the unnormalized noise-sum threshold \(k=2|g|\) is:

| Mass | kg | \(k=Gm^2s^2/(\hbar d^3)\), s⁻¹ |
|---|---:|---:|
| 10 fg | \(10^{-17}\) | \(6.32891937\times10^{-10}\) |
| 10 pg | \(10^{-14}\) | \(6.32891937\times10^{-4}\) |
| 1 ng | \(10^{-12}\) | \(6.32891937\) |

The normalized Ramsey decay-sum threshold is \(r_A+r_B\ge4|g|=2k\).
This factor of two is a normalization distinction, not a physics discrepancy.

The PDF of arXiv:2603.26075v2 **page 8, equation (36)** prints approximately
6 Hz for **10 fg** with these distances. Direct SI substitution gives the first
row instead: the printed coefficient is about \(9.48029177\times10^9\) times
too large. A 1 ng mass gives a rate near 6 in this formula. Neither Hz versus
angular frequency nor the 1/99 geometry correction explains the discrepancy.
This is an apparent numerical/unit error in the displayed benchmark; it does
not disprove the symbolic inequality. No author correction is asserted.

**Feasibility is not certified.** The smallest displayed benchmark has an
inverse rate of approximately 50 years. That is an information-scale warning,
not a proved integration time or sample-size requirement. The stronger masses
also require controlled coherent preparations and suppressed cross-talk.
Recent atom-interferometry work remains a proposed platform
([Howl et al., PRA 114, 023306](https://journals.aps.org/pra/abstract/10.1103/l62d-gz5c)).
No observed gravitational entanglement or qualifying experimental lock is
claimed here.

The same PDF's page 21, equation (129), writes \(e^{if t}\) although (127)
is \(\dot\rho=f\rho\) and (128) defines a complex coefficient. We implement
\(e^{ft}\), verified directly against the generator. The printed extra i would
turn dephasing into an incorrect oscillation. This is a notation error, not a
reason to abandon the correctly derived generator.

## 5. A memory obstruction to over-broad interpretation

Time-local differential equations do not alone imply CP-divisible dynamics.
This distinction is established in open-system theory
([Breuer–Laine–Piilo](https://arxiv.org/pdf/0908.0238), equations (6)–(9)).
Here is a controlled conventional boundary example:

\[
H_{SE}=Z_S\otimes\operatorname{diag}(-1,0,1),\quad
\sigma_E=\operatorname{diag}(1/4,1/2,1/4),\quad
q(t)=\cos^2t.
\]

Tracing the three-state environment gives the exact qubit channel
\(\rho_{01}(t)=q(t)\rho_{01}(0)\), with populations unchanged. Every full
map from preparation time is CPTP. On intervals where \(q>0\), its time-local
equation is \(\dot\rho=\tan t\,(Z\rho Z-\rho)\). The coefficient is negative
on \((\pi/2,\pi)\). Between \(2\pi/3\) and \(5\pi/6\), the intermediate
coherence multiplier is 3; its normalized Choi matrix has eigenvalue −1.
Thus the intermediate map is not CP although the complete histories are physical.
The environment can equivalently be a retained classical random phase label.

This precisely defeats the inference “time-local alone implies a positive-rate
Lindblad generator.” It does **not** construct a viable non-Markovian classical
gravity model or refute the discriminator inside its positive-rate class. A
finite experiment can bound memory on its tested repertoire; it cannot promote
that bound into a universal Markovian premise.

## 6. What a gravitational entanglement signal would and would not decide

[Di Biagio](https://arxiv.org/html/2511.02683), accepted PRD 1 September 2026,
distinguishes circuit mediation/state-space factorization from spacetime causality.
[Gundhi–Infantino–Bassi](https://arxiv.org/html/2604.19696v3) identify matter
diffusion and the chosen mode partition as confounders in the specific
Aziz–Howl construction. These are different arguments, not competing experimental
detections. This audit does not adjudicate every gravitational-field ontology.

The safe conclusion of the operational certificate is **exclusion of the
specified non-entangling positive-rate effective dynamics**, with the stated
interface and error bounds. It cannot certify a unique theory of reality,
quantum ontology of every gravitational degree of freedom, a generated observer,
consciousness or the cosmological origin of physical laws.

## 7. Reproduction and verdict

Run from the repository root:

```bash
bash research/recon01/reproduce.sh
```

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. Results, exact rational determinant
and interval checks, matrix-unit identities, independent matrix-exponential
comparison, finite trajectories, the memory dilation and SI estimates are in
[results.json](results.json). [MANIFEST.json](MANIFEST.json) pins the files.
Numerical residuals support implementation correctness; the algebraic theorem
and exact geometry identity carry the whole-family statements. Floating-point
tolerances are diagnostics, not certified physical uncertainty enclosures.

Execution returned `PUBLISHED_DISCRIMINATOR_AUDIT_PASS_NO_NEW_GRUT_LAW`:
8 disclosed fixtures, 128 generator matrix-unit cells and 896 time/basis cells.
Both the partial-transpose generator identity and kernel projection had zero
matrix residual. The independent exponential comparison had zero residual at
the reported precision. The finite-step slope diagnostic differed by at most
\(3.751\times10^{-7}\,\mathrm{s}^{-1}\). The smallest Schur eigenvalue was
\(-5.67\times10^{-16}\), within floating-point tolerance. Four synthetic interval
cases returned exclusion, inconclusive, inconsistent calibration and invalid
input respectively. No experimental result is inferred from these checks.

**Verdict: RETAIN THE DISCRIMINATOR; NO NEW GRUT LAW ESTABLISHED.**
The useful output is a defined experiment/comparator calculation with exposed
memory and calibration boundaries. The searched constructions have not supplied
an admissible new joint-generating law; adopting their published equations would
spend a card without establishing that novelty. Card 1 remains author KILL,
pending independent review; two slots remain. This completed packet does not
reopen a foundations queue or call itself an external review.

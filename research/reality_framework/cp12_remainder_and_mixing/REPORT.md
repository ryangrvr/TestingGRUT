# CP-12: graph mixing excludes a ground phase with vanishing degree defects

**Result:** at CP-12's frozen couplings and \(Q=2N\), the complete
zero-temperature ground state cannot concentrate on four-regular graphs,
or on graphs with only \(o(N)\) total absolute degree defects. This is a
full-Hamiltonian variational result, including graph mixing and all
integer electric fields. It does not use an unbounded cluster remainder.

This sharpens the question posed by the even-width review at
[`e7a6535`](https://github.com/ryangrvr/TestingGRUT/commit/e7a6535e849d876a0f7b4611c0a2d98ae6437148).
Instead of trying to certify a difference near \(10^{-9}\) between two
regular graph blocks, we can exclude their common nearly regular sector
as a dominant ground phase. The same conclusion includes regular
nonspatial competitors. It does not determine which alternative geometry
the true ground state selects.

The full proof is in [DERIVATION.md](DERIVATION.md). This directory is an
additive follow-up; frozen CP-12, earlier candidates, canonical GRUT-RAI,
PR #12, and Stage-3 Card 3 are not changed.

## The physical mechanism behind the result

CP-12's neutral reconnection operator is a positive connection Laplacian.
For a fixed matter/field configuration, its graph fiber is exactly a
Johnson graph. A uniform superposition over compatible connections has
zero rewiring energy. Holding geometry fixed retains the Laplacian's
positive diagonal while suppressing that coherent cancellation.

Two bounds make this distinction quantitative. Gauss law and the fermion
anticommutator give every complete fixed-graph block at \(Q=2N\) the
lower bound \(E_G/N\ge\ell_N\), with

\[
\ell_\infty=\frac{353-\sqrt{128321}}{160}
=-0.03262004144948\ldots .
\]

A physical trial state coherently combines local vertex-edge charge
matchings and all compatible graph configurations. It satisfies both
Gauss laws and lies in the rewiring kernel. Its Rayleigh energy density
approaches

\[
u_\infty=49/32-\sqrt{21/8}
=-0.08893517460197\ldots .
\]

The difference is about **0.056315 per vertex** in the frozen model's
units. This is an extensive energy advantage, not a further tiny
cluster coefficient. The trial state is used only to bound the ground
energy; it is not assumed to be a prepared universe or the actual ground
state.

Regular-to-nonregular coupling is much smaller than the full off-diagonal
norm. One edge swap has at most five regular predecessors at degree four:

\[
\|P_4H(I-P_4)\|\le\frac\Gamma M\sqrt{5Q(M-Q)}=O(N^{-1/2}).
\]

Combining this with the extensive energy separation gives

\[
\boxed{\Pr_{\rm ground}(\text{every degree is four})=O(N^{-3}).}
\]

Exact regularity is an overly strict physical target. The stronger
weighted-rewiring estimate uses
\(D(G)=\sum_i|d_i(G)-4|\). It gives

\[
\boxed{\Pr_{\rm ground}[D\le K_N]\to0\quad\text{for every }K_N=o(N).}
\]

One conservative quantitative corollary is
\(\liminf\langle D\rangle/N\ge39/400000\). This is a deliberately weak
model bound on a positive defect density, not an observed property of
space. The proof and exact rational comparisons are provided.

## Executable controls and finite certificates

The new trial compression and local bound were checked against complete
physical matrices, rather than inferred from cluster coefficients:

| Control | Size | Result |
| --- | ---: | --- |
| Coherent matching matrix and Gauss constraints | 30, 872, 3,975 states | Predicted matrix error at most \(4.30\times10^{-14}\) |
| Degree-four central-star lower bound | Complete 126-state graph block | Local minimum agrees within \(5\times10^{-15}\) |
| Star sum and electric/removable diagonal split | Three complete graph blocks | Reconstruction error at most \(4.26\times10^{-15}\) |
| Nearly regular weighted and crossing bounds | All 455 graphs at \(N=6,Q=12\) | Pass |
| Matching-count ratio used for large certificates | Exact counts through \(N=128\) | Pass |

Finite probability bounds for **exact** four-regular ground weight use
80-digit outward arithmetic, not floating-point eigensolver precision:

| \(N\) | Upper bound on exact four-regular ground probability |
| ---: | ---: |
| 32 | 0.010404 |
| 64 | \(4.0649\times10^{-4}\) |
| 256 | \(4.3111\times10^{-6}\) |
| 1,024 | \(6.0908\times10^{-8}\) |
| 2,048 | \(7.4818\times10^{-9}\) |

The conservative nearly regular certificates cover all graphs with
\(D\le K\), not only any selected lattice examples:

| \(N\) | \(K\) | Upper bound on ground probability of \(D\le K\) |
| ---: | ---: | ---: |
| 65,536 | 16 | 0.18405 |
| 262,144 | 22 | 0.021251 |
| 1,048,576 | 32 | 0.004561 |
| 16,777,216 | 64 | 0.000337 |

These rows are analytical certificates evaluated at large \(N\). They are
**not** diagonalizations of million-event Hamiltonians. The files retain
the explicit trial energy, sector lower bound, gap, and coupling estimate.

## Complete small systems show why graph-block rankings are insufficient

The finite full-law controls use \(Q=N\), a separate bond density from the
four-regular theorem. They include every labelled graph and all regulated
physical matter/field states. At \(N=4,Q=4,L=3\), the best single graph
block is a regular square, but the coherently mixed nonregular sector has
lower energy:

| Energy, including the actual rewiring diagonal | Value |
| --- | ---: |
| Best regular graph block | −0.5116660370 |
| Best nonregular graph block | −0.4803773622 |
| Full nonregular compressed sector, retaining its graph mixing | −0.5631298716 |
| Complete Hamiltonian ground | −0.5985522825 |

The full ground weight on regular graphs is 0.28956. At
\(N=5,Q=5,L=3\), all 427,986 states are retained, and that weight is
0.09487. Cutoffs 2 and 3 give full-ground energy differences below
\(1.9\times10^{-10}\). These numerical comparisons support the mechanism;
they are not the proof of the \(Q=2N\) large-system result.

There is also a positive dynamical result: for any coherent initial state
entirely in an exact regular sector, and norm-one graph-diagonal source
and detector operators, the retarded-response error against compressed
regular evolution is at most \(4b_N^2T^2=O(T^2/N)\). Complete finite-law
response tests pass. Prepared regular geometry can be stable on fixed
time intervals even while its ground-state selection fails.

## Remainder work and the incoming execution audit

The exact committed even-width script fails before execution with a
`SyntaxError` in its output dictionary. `even_width_repaired.py` quotes
the dictionary keys and changes two mistaken assignments to dictionary
entries; its embedding and energy calculations are unchanged. The repaired
script executes and reproduces both eight-edge increments. The original
pinned script is not represented as having passed. The account and minimal
diff are stored beside the numerical output.

Complete periodic \(C_4,C_6,C_8\) graph blocks independently reproduce
their finite linked sums at cutoffs 3 and 4. This validates finite-host
tails, **not** the infinite square-lattice tail. A separate complete
1,032-state Schur-resolvent control reproduces graph mixing with a
controlled Chebyshev polynomial remainder and a retained floating-point
error floor. It is a finite-cutoff validation, not an interval-certified
large-system eigenvalue calculation. These controls remain useful, but
the phase exclusion above bypasses their unresolved infinite-host tails.

## Reproduction

Run from the repository root with Python, NumPy, SciPy, and NetworkX.
The parent commit supplies the frozen model and cluster data. Set BLAS
threads to one for predictable resource use. Every recorded result here
was executed; the fast verifier also checks their mutual identities.

```bash
python research/reality_framework/cp12_remainder_and_mixing/verify.py
python research/reality_framework/cp12_remainder_and_mixing/star_bound_check.py
python research/reality_framework/cp12_remainder_and_mixing/variational_exclusion.py
python research/reality_framework/cp12_remainder_and_mixing/near_regular_exclusion.py
python research/reality_framework/cp12_remainder_and_mixing/even_width_repaired.py
python research/reality_framework/cp12_remainder_and_mixing/regular_graph_bound.py
python research/reality_framework/cp12_remainder_and_mixing/mixing_fibers.py --n 4 --q 3 --L 3
python research/reality_framework/cp12_remainder_and_mixing/schur_remainder.py
python research/reality_framework/cp12_remainder_and_mixing/finite_host_tail.py
python research/reality_framework/cp12_remainder_and_mixing/regular_full.py --n 4 --d 2 --L 3
python research/reality_framework/cp12_remainder_and_mixing/regular_full.py --n 5 --d 2 --L 3 --no-response
```

The last command is the larger many-body calculation. The verification
manifest records dependencies, hashes, sizes, and the executed parameter
sets. Floating-point residuals support matrix controls; outward arithmetic
and inequalities support the analytical certificates. They are separate
forms of evidence.

## What this means for the research

We have answered a full-Hamiltonian selection question without choosing a
new action or inferring a phase from a sign-changing expansion. Frozen
CP-12 does not select a ground ensemble approaching exact degree-four
organization with vanishing absolute degree defects.

A spatial phase with a **nonzero** defect density remains possible. So do
other bond densities and metastable prepared geometries. Degree
irregularity alone is not nonspatiality; a manifoldlike sprinkling need
not be regular. The next unresolved physical comparison must therefore
allow irregular spatial and nonspatial competitors in the same complete
state prescription, and compare their propagation and macroscopic growth.

The microscopic challenge is now more concrete: graph coherence has an
extensive energetic benefit, while conditional matter propagation can
remain stable on prepared geometry. A successful construction must make
those two properties coexist in a physically selected spatial phase.
The present theorem supplies a benchmark that such a construction must
survive. It supplies neither a replacement law nor a completed theory of
reality.

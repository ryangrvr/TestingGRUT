# Sources and attribution limits

## Inspected project material

* S4 packet and executable model at
  https://github.com/ryangrvr/TestingGRUT/commit/c6ac0c68480d4164ad620bbf2a4b3d8c65ac85b1 .
  S5 imports its immutable Python complete-model construction for
  independent numerical comparisons, and uses its conditional
  protection theorem for the analytic joined protocol.
* Published S5 charter:
  https://github.com/ryangrvr/TestingGRUT/commit/64378573e8251224dd563928667e961a0effb3da .
* Published P1 charter:
  https://github.com/ryangrvr/TestingGRUT/commit/364781b841121e3169041b121e661048fd793524 .
* The user's supplied S4 assessment and evidence ladder. Its referenced
  reviews/grut_review_report.md was not available here and is not
  claimed as independently inspected.

## Primary methodological source

Andrew M. Childs, Yuan Su, Minh C. Tran, Nathan Wiebe and Shuchen Zhu,
"A Theory of Trotter Error", arXiv:1912.08854v3; published as
"Theory of Trotter Error with Commutator Scaling", Physical Review X
11, 011020 (2021): https://arxiv.org/html/1912.08854v3 . Section 5.1,
Eq. (150) was checked directly for the two-term unitary Strang bound
and its 1/12 and 1/24 prefactors. The project-specific commutator
constants are derived here, not imported as a theorem about GRUT.

Standard tools include bosonic second quantization, discrete sine
diagonalization, natural orbitals, two-mode Fock transformations,
unitary error telescoping, Duhamel and projector triangle inequalities.
No general algorithmic novelty, exhaustive priority review, physical
constant selection or distinctive empirical prediction is claimed from
using these tools. The literature searches did not supply a uniform
formation theorem for this history-retaining model; none is inferred.

All project-specific formulas used in acceptance are stated in
MATHEMATICS.md and PROBABILITY_BOUND.md. Computational comparisons by
the same Builder do not constitute independent external validation.

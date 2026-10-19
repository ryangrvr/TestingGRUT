# C2 response and corrections before the B control

19 October 2026. Follow-up to `0b6df6f090f7d02d3de2cff132e30ce53959d162`.
This records the owner's withdrawals and answers the question about P1. It
does not execute B Step 1, select a new state, or advance to Steps 2–3.

**Answer:** C2 identifies the appropriate observer-state object, but P1 does
not currently add a state-selection condition. There is no stronger, hidden
closure condition in the `0b6df6f` derivations. That commit formulates the
selection problem; it does not solve it. The present project outcome remains
the map of established results and limitations (option C).

## What the archived P1 actually says

At the cited commit, [GRUT III P1](../../../GRUT_III/01_PRINCIPLES.md)
is labelled an identification, with zero additional postulate content. It
starts with a global state and assigns an observer the triple

\[
\bigl(\mathfrak A(O),\;\omega|_{\mathfrak A(O)},\;\sigma_t^\omega\bigr).
\]

The bath is the commutant; identification with the algebra of the causal
complement additionally requires Haag duality. P1 expressly supplies no
influence functional or dynamics. [NEXT_QUESTION.md](NEXT_QUESTION.md) asks
what extra condition could select a state or cycle; it supplies none.

The restriction map

\[
\omega_O(A)=\omega(A),\qquad A\in\mathfrak A(O),
\]

exists for every state on the encompassing algebra. Consequently the
instruction to restrict a state does not choose that state. In the free-field
setting, the vacuum and smooth coherent excitations both admit this
restriction and the same commutant assignment, although their local field
expectations differ. Thus P1 alone cannot fix even the local one-point
function. This is an elementary limitation of the definition, not a new
quantum-cosmology calculation.

Likewise, for a faithful normal state the KMS property with respect to its
own modular flow is automatic. It does not select the state or set a
temperature in physical time. Equating that flow to a separately specified
geometric observer evolution would add content; P1 has no such general
condition for an inflating subregion. Its special geometric examples are
already labelled imports in the archive.

There is also a distinction between the two readings of C2. Treating the
observer restriction as the object of prediction follows P1. Prohibiting any
global state and defining a state only on observer algebras would be a new
construction choice: P1 itself uses a global state. That choice still needs
a rule that assigns the local state, and compatibility rules for overlapping
observer algebras. It does not by itself choose one of the proposed lapse
cycles.

## What the density-matrix paper already contains

The full reference is **Victor Ivo, Yue-Zhou Li and Juan Maldacena, “The no
boundary density matrix,” JHEP 02 (2025) 124**, published 19 February 2025,
DOI [10.1007/JHEP02(2025)124](https://doi.org/10.1007/JHEP02%282025%29124),
[arXiv:2409.14218v2](https://arxiv.org/html/2409.14218v2).

The introduction, eqs. (2)–(4), starts from the full Hartle–Hawking functional
and integrates over the unobserved configuration. At the trace saddle it
matches the fields and their conjugate momenta on the two sheets in that
region. These matching conditions follow from tracing and stationarity;
they are already part of the published prescription. They are not an extra
GRUT condition obtained by renaming the complement a bath. The paper also
explicitly leaves normalization/selection issues open in its introductory
discussion. Restriction alone does not repair them.

Strictly, the continuum local algebra used in P1 has no intrinsic partial
trace or canonical tensor factorization into inside and outside Hilbert
spaces. Its state is a positive normalized functional. A literal subregion
density matrix needs an explicit regulator or suitable split construction;
gravity also requires a prescription for the subregion observables and
constraints. This qualification does not invalidate the paper's
semiclassical path-integral kernels. It prevents silently replacing P1 by a
finite-dimensional tensor-product model. See Witten,
[arXiv:1803.04993v6](https://arxiv.org/html/1803.04993v6), sections 6.4–6.5
and 7 (with the formal path-integral example in section 5.2).

An extension of C2 must therefore name a condition on the local state or
gluing prescription that is **not** already restriction, trace matching,
positivity, normalization, or modular KMS. No such extension has been derived
here. Reproducing the paper may be a useful later control, but does not
itself create a new principle.

## Full Obadia reference and the withdrawn D2 claim

**Nathaniel Obadia, “How hot are expanding universes?”, Physical Review D
78, 083532 (2008)**, published 23 October 2008. DOI
[10.1103/PhysRevD.78.083532](https://doi.org/10.1103/PhysRevD.78.083532).
Preprint [arXiv:0804.2890](https://arxiv.org/abs/0804.2890).

For the flat-FLRW conformal massless setup, eqs. (20)–(22) give the relevant
short-separation coefficient. For a comoving observer it yields

\[
T_{\rm coefficient}=\frac{\sqrt{H^2+2\dot H}}{2\pi}
=\frac{H}{2\pi}\,[1-\epsilon+O(\epsilon^2)],
\qquad \epsilon=-\dot H/H^2.
\]

The coefficient alone is not a thermal-state theorem: the paper separately
requires the stationary thermal Wightman structure. This is the prior source
for the first-order local coefficient; it is not a claim that every result
in v11/v14, every field, or the closed-slicing problem is covered by that
formula. D2 is withdrawn as a claimed new result and retained only as an
appropriately scoped consistency check. No novelty certificate is issued.

## Source corrections needed before registering Step 1

The owner's proposed ordering is retained: the published-prescription
control comes first. The source audit finds the following distinctions.

| Proposed case | Correct scope for the control |
|---|---|
| (a) FLT, 1703.02076 / 1705.00192 | Fix the Dirichlet data, positive-real lapse prescription and regular perturbation solution. The quoted inverse-Gaussian result belongs to that specified prescription. |
| (b) Diaz Dorronsoro et al., 1705.05340 | The contour covers the whole real axis and bypasses zero below; section III gives a single inhomogeneous scalar calculation, eqs. (3.20)–(3.23). It is not the positive half-line. |
| (b), companion 1804.01102 | The main exact model is nonlinear biaxial Bianchi IX with a mixed initial size/momentum prescription and a closed lapse contour. It is a separate model check, not directly a third scalar width. |
| (c) Di Tucci–Lehners, 1903.06757 | A specified family of Robin conditions gives stable HH saddles. Off-shell initial geometries need not have zero size. Arbitrary Robin and Neumann data are not synonymous with this family. |
| Yamada, 2511.09621 | Separate regulator/initial-state prescription: integrate with a finite-width initial-size Gaussian before taking its zero-width limit; include an initial perturbation wavefunction. The positive lapse version gives a tunneling wavefunction with suppressed perturbations, while a whole-line contour below zero also recovers HH. It does not establish “stable HH only.” |

The Yamada journal reference is **Masaki Yamada, “New regularization scheme
for the wave function of the Universe in the Lorentzian path integral,”
Physical Review D 113, 063547 (2026)**, published 18 March 2026,
DOI [10.1103/z59k-v3l8](https://doi.org/10.1103/z59k-v3l8). The arXiv v2
text, checked here, was revised 24 August 2026. Claims about its prescription
must specify the version and order of limits.

Accordingly, a three-prescription script can compare the specified FLT,
whole-line HH and Di Tucci–Lehners Robin controls. Yamada should be marked a
distinct prescription, not used as an interchangeable citation for the
third result. These are source identifications, not reproduced widths or
an adjudication of the full contour dispute.

Before computation, fix the scalar action, harmonic normalization and label,
final size, boundary terms, cycles, regulator limits and amplitude versus
probability convention. Exact de Sitter admits a scalar spectator control;
it should not be silently called an inflaton curvature perturbation. If

\[
\Psi_\alpha(z_1)\propto
\exp[-(B_\alpha+iC_\alpha)z_1^2/(2\hbar)],
\]

then suppression means \(B_\alpha>0\); its normalized Gaussian probability
has variance \(\hbar/(2B_\alpha)\). For \(B_\alpha<0\), report an
inverse-Gaussian coefficient, not a physical variance. A conjugate-saddle
sum can have an oscillatory interference factor; distinguish its Gaussian
envelope from its complete wavefunction. A sign match alone does not count
as reproducing a contour integral, and the expected sign must not be put
into the numerical boundary conditions as the result being tested.

Finally, C1 is not logically identical to (a). “No boundary data of any kind”
does not determine a lapse integration cycle; Dirichlet zero size is itself
boundary data. If C1 is **defined** to be (a), its inverse-Gaussian failure
is a conditional termination of that definition, not of every meaning of
closure. Similarly, normalizability is an existing admissibility condition,
not a new GRUT selection principle.

## Decision and execution boundary

- D2 novelty, a sign-only curvature falsifier, and the claim that all
  sequestering formulations force positive spatial curvature remain
  withdrawn as requested.
- C2 as presently stated adds an observer interpretation; it does not yet
  supply the missing selection condition. There is no alternate derived
  GRUT condition to substitute for it.
- Step 1 has not been executed in this follow-up. No three widths, contour
  reproduction, density-matrix spectrum, KSW bound or new observable is
  claimed. Steps 2–3 remain conditional on the prescribed sequence.
- This bibliographic/definition clarification stays on the existing branch;
  it opens no new research branch. The one-branch-per-step and at-most-two
  branches/day limits are retained for executing B.

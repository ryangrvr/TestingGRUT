# Post-campaign analytical extension: state-weighted positive residual

This derivation follows the stopped S7 campaign. **It was not used in any
S7 trajectory or error ledger.** Its small independent matrix controls are
additional tests, not retroactive permission to pass a frozen gate or a
claim that six-event formation is solved.

## Identity before any time approximation

Use an exact common receiver isometry W and its complete word/Fock lift V.
Let P=VV*, C=V*BV and D=(1-P)BV. For a word of k non-root attachments set

\[
R|_{u}=\sqrt{M-k}\,I,\quad M=N-1,\quad
\sigma=\|(1-WW^*)c\|.
\]

Unique parent words and the outside-photon vacuum give the operator identity

\[
\boxed{D^*D={g^2\sigma^2\over N^2}R^2.}
\]

To see this, annihilation stays inside P. A child's outside component can
arise only by append from its unique parent. Different child words are
orthogonal. Outside creation acts as an isometry times sigma on every
inside Fock vector, since the outside mode is empty. Summing remaining
type multiplicities for a k-word gives M-k. This proves the identity for
arbitrary coherent word packets, not only diagonal geometry mixtures.

Therefore `||Dx||=(g*sigma/N)*||Rx||`. The S7 operator estimate replaces
`||Rx||` by `sqrt(M)*||x||`; it can be loose late in formation.

## A step bound must pay changes during that step

It would be incorrect to multiply the initial weighted residual by h and
ignore the subsequent coefficient evolution. C changes k by one. On each
edge, the multiplier in `[R,C]` has magnitude
`|sqrt(r)-sqrt(r-1)|<=1`. Entrywise absolute domination by |C| and the same
complete bosonic emission norm used in S5/S7 give

\[
\|[R,C]\|\le G,\qquad G=g\|W^*c\|.
\]

One-particle orbital rephasing makes the coefficient contact components
nonnegative, so the complete absolute-entry emission norm obeys this G
bound as well. No geometry transitions or reverse edges are removed.

The unitary commutator identity then yields

\[
\|R e^{-iCs}x\|\le\|Rx\|+|s|G\|x\|.
\]

Substitute the exact D identity into the physical-H Duhamel integral:

\[
\boxed{
\|e^{-iBh}Vx-Ve^{-iCh}x\|
\le {g\sigma\over N}
\left(|h|\|Rx\|+{h^2G\over2}\|x\|\right).
}
\]

The minimum of this bound and S7's operator bound is valid. Even a fully
assembled initial coefficient vector (Rx=0) retains an h^2 allowance:
reverse detachment can make append possible during the step. Setting its
whole-step residual to zero would discard precisely those returns.

The weighted norm is a positive sum
`||Rx||^2=sum_u (M-|u|)*||x_u||^2`; it needs no subtraction of total and
captured action norms. A future implementation must enclose this sum, the
contact defect, coefficients, polar-frame interpretation and every other
error exactly as in S7. The displayed identity assumes exact isometries;
it does not excuse those numerical costs.

## Verification and scope

`weighted_control.py` builds full small physical receiver matrices using
the independent immutable S5 Python operator and Fock lifts, tests D*D,
the commutator bound and full emission
wavefunctions for random and fully assembled packets. It uses ordinary
double precision for these diagnostics. It does not execute a new primary
formation family or apply the proposed bound retrospectively to S7.

This gives a concrete next residual construction, not a rank-growth or
formation theorem. The h^2 term, projection errors and growing Fock space
remain. A next numerical campaign must freeze its implementation and gate
before rerunning recovery and attempting T=144.

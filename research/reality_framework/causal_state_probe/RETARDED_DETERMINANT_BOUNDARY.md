# Retarded determinant does not weight causal histories

Date: 2026-10-11. Additive analytical and executable control. No CP-07, proposed fundamental law, altered frozen packet, PR modification or continuum-emergence claim.

## Research finding

**Exact finite-order obstruction.** Any retarded matrix R on an acyclic finite causal order has a topological labelling in which R is strictly triangular, hence nilpotent. Therefore for arbitrary complex z,

\[
\det(I+zR)=1,\qquad {\rm tr}\,R^k=0\quad(k\ge 1).
\]

More generally, for any nonzero constant a, \(\det(aI+zR)=a^N\). This is independent of order statistics, links, height, and the commutator spectrum. It holds for *any* retarded weighting with zero diagonal that vanishes on unrelated/advanced pairs, not just the link kernel.

It defeats a tempting candidate for matter-driven selection: taking a relative history weight from \(\det(I+zG_R)^{-1/2}\). This weight cannot distinguish *any* two finite causal histories under the stated prescription. In addition, a nonsymmetric retarded operator is not by itself a legitimate Hessian of a real scalar Euclidean quadratic action; using its determinant as a path-integral weight requires a separately justified microscopic action and measure.

**Counter-control: supplied symmetric Gaussian action.** If one instead defines \(S=(R+R^T)/2\) and \(K_E=I+zS\), for sufficiently small real z this is a real positive-definite Hessian and its Gaussian integral yields \((\det K_E)^{-1/2}\). It has nontrivial history dependence:

\[
\log\det(I+zS)=-(z^2/2){\rm tr}\,S^2+O(z^3)
=-(z^2/4)\,L_{\rm links}+O(z^3).
\]

This is a *different, additionally supplied* action: symmetry/positivity do not derive it from causal order, nor give it the Lorentzian gravitational and matter correspondence required for GRUT. The choice of kernel, coupling z, measure, counterterms, and ensemble are independent inputs.

## Numerical control (NumPy, 180 events)

Reproduced the reviewed 4D sprinkling, and generated the same 90-by-90 random two-layer control with identical relation count, same seeds and identical bare link prescription. With z=0.01 for the symmetric Gaussian:

| Quantity | 4D sprinkled interval | Height-2 random order |
|---|---:|---:|
| Event count | 180 | 180 |
| Relation count | 1413 | 1413 |
| Link count | 897 | 1413 |
| log det(I+z R) for z=0.013, -0.08, 0.4 | 0 (numerical tolerance) | 0 (numerical tolerance) |
| log det K_E | -0.0224455689292 | -0.0353575261699 |
| Gaussian log-weight, -1/2 log det K_E | +0.0112227844646 | +0.0176787630849 |
| Minimum eigenvalue of K_E | 0.937668125 | 0.917469227 |

For this *illustrative* symmetric kernel, a single random layered impostor receives a larger Gaussian factor than the single sprinkled geometry. The ratio is only about 1.00648; it says nothing about total dominance after summing complete history classes with a common measure. It does demonstrate that merely replacing a retarded matrix with a symmetric Gaussian kinetic operator does not, by itself, earn manifold selection.

Independent positive control: a directed three-cycle (not a causal order) has \(\det(I+zR)=1+z^3=1.027\) at z=0.3. The determinant routine is not trivially hardwired to output one.

## Scientific boundary and next step

- The original rank-two control was a property of complete bipartite links. The independent same-relation-count random layered control achieves commutator rank 180, so rank is not a geometry selector.
- The determinant theorem is stronger: **a strictly retarded determinant offers no causal-history discrimination at all**, even when the commutator rank varies.
- A legitimate matter-induced history amplitude requires a fixed *physical quadratic or interacting action*, its state and boundary prescription, calibrated density/normalization, treatment of zero modes, regulator and measure, and a correct probability/amplitude comparison between *classes* of orders—not only two selected orders.
- A quantum matter state on one history (including an SJ state) is not automatically a relative amplitude of that history in a sum-over-histories.

The next useful constructive test is thus **not** another rank, pair-count or determinant check. It is to derive a physical kinetic/matter operator from a proposed microscopic rule, then test whether its induced relative weights suppress nonmanifoldlike histories without assuming the desired spacetime or its action. If the rule supplies a preferred metric, signature, dimension, or geometric ensemble, list those as supplied premises. The present calculation does not choose that rule.

## Reproduction and provenance

Run from this folder:

\`\`\`bash
python retarded_history_weight_no_go.py
\`\`\`

Expected final line: \`PASS: retarded determinant geometry-blind; symmetric action responds but is supplied\`.

Parent reviewed state and source scripts: [REVIEW.md](REVIEW.md), [supplied_probe.py](supplied_probe.py), [adversarial_rank.py](adversarial_rank.py). Added only in this additive research directory; CP-06 and older frozen work are untouched.

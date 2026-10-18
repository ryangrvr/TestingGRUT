# Independent review route

No independent reviewer has validated this packet. One assistant produced
the construction, proof, challenge and controls with shared code ancestry.

## Load-bearing steps to try to break

1. **Exact half-line sector.** Reconstruct h from B1's local operators.
   In particular the root co-pair energy is 2-A, bulk co-pairs have energy
   2, fuel energy Omega=2-2A, and every relative separation is retained.
2. **Weighted graph.** Check the rectangular weighted Schur norm, completeness,
   the contraction map and its domain with unbounded D. Verify that the
   invariant graph equals the physical spectral subspace, using its
   orthogonal-complement form estimate. A bounded graph generator must
   not be called Hermitian in the original P inner product.
3. **Growth, not phase approximation.** Reproduce the parity cancellation and
   bound on the EXACT correction E. Check logarithmic growth of the weighted
   exact generator, graph reconstruction and the inverse metric on the
   ORIGINAL fuel vector. No approximation may be exponentiated in place of it.
4. **Joint embedding and defect.** Check validity for separated strings,
   the inherited cut diagonal, and the incoming AND outgoing boundary
   drives in (13). The comparator is on the same full physical space and
   is unnormalized. The physical U_L is unchanged throughout Duhamel.
5. **Initial high weight.** Verify the one-source gap via the Schur resolvent,
   variance J^2, and input error 2sqrt(w). This component must not be charged
   again at every time in the integral.
6. **Original event.** Check R_L Phi's amplitude is a_L, not sqrt(a_L), since
   it is a product of TWO identical spatial packets; the unknown qubits
   need not be a product state. Color factors as I_4, so the estimate and
   exact successful-output channel hold with an external entangled reference.
7. **Remaining source inequality.** Neither PROOF.md (18) nor SOURCE_ESTIMATE.md
   (22) is proved. They must remain open even if every preceding step passes.
8. **Short launch inheritance.** B1's finite-path bounds are uniform in path
   length; pass to the half-line using strong-resolvent / finite-support
   convergence on the same closed quadratic-form domain. The short-time
   bound is independent of its full-geometry extension. Do not extrapolate
   it beyond theta=4 or reinterpret it as distant output.

The shortest complete mathematical route is PROOF.md sections 2-5, then
SOURCE_ESTIMATE.md. CHALLENGE.md lists scope and inference failures.

## Reproduction and chronology

* Parent: `9256c11ba554a0d080b695dec49fc413866f6455`.
* Charter/pins freeze: `abf838071590232053581a2ef3fa2098a335d90a`.
* Proof/source freeze: `f8f5dfcf0a24f42e0267dcfdf7adc98475d200f9`.
* Then one control run and one identical-source replay.
* SOURCE_ESTIMATE.md is a later analytic translation to unfiltered source
  observables, not an altered physical gate or a new primary run.

From repository root, integrity is checked by

```sh
python research/reality_framework/grut_b2_source_output_20261018/verify.py
```

The frozen controls can be reproduced separately with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python research/reality_framework/grut_b2_source_output_20261018/controls.py /tmp/b2_so_review
```

They check rational constants and matrix embeddings only. They do not
evolve u(t), evaluate source launch or the operating probabilities, compute
the infinite graph numerically, or validate the all-size mathematics.
Matrix arithmetic is ordinary double precision without interval enclosures.

## Relation to known methods

Exact block-diagonalization and low-energy analysis are established
methods; see Bravyi, DiVincenzo and Loss, *Schrieffer-Wolff transformation
for quantum many-body systems*, Ann. Phys. 326 (2011), 2793-2826,
https://arxiv.org/abs/1105.0675. We do not import its many-body theorems as
a proof for this unbounded half-line space. Our fixed-point, weighted
constants, reconstruction and full-H joint defect are explicitly derived.

The user's comparison paper, Blanco de Paz and Huidobro, *Bound states in
the continuum in subwavelength emitter arrays*, Phys. Rev. Research 5,
033108 (2023), https://doi.org/10.1103/PhysRevResearch.5.033108, concerns
emitter arrays. Its existence of nonradiative states is context, not a
diagnosis or proof about the present fuel spectral measure. No external
paper supplies the unresolved source lower bound here.

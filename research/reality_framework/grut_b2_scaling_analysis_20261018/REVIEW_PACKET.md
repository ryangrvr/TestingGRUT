# Short review route

No independent mathematical reviewer has inspected or approved this packet.
All derivations and challenges are from the same assistant. The code reuses
B1/B2's scalar catalog and colored builder. Agreement is not independence.

## What to try to falsify

1. **The all-L band construction.** In DERIVATION.md (1)-(8), check Q's true
   minimum diagonal, the resource-count hopping norm, the common circle,
   the contour estimate for Pi-P, and the exact canonical isometry. The
   configuration-distance claim must include all intermediate charge/fuel
   states and coherent source/receiver sets. It is not a spatial lower bound.
2. **The preparation-weight lemma.** Check min-max's high edge, exact initial
   energy variance, and pure-state trace distance. In particular do not
   confuse virtual coordinate occupation with high spectral weight.
3. **The bulk velocity bounds.** Reconstruct the simultaneous-translation
   fiber, its phase gauge, the positive isolated even ground state, and
   Hellmann-Feynman derivative. Check both u0's lower bound and u1's upper
   bound used in the post-freeze estimate (18). No finite relative cutoff
   or effective hopping approximation is allowed in those proofs.
4. **The information theorem.** Test color order on the extended path with
   endpoint fuel sites, including states with both a fuel and plus carrier
   at the endpoint event. The theorem is about arbitrary two-qubit inputs
   and an outside reference, but only on a fixed path with endpoint fuels.
5. **The proposed route's missing inputs.** Equations (19)-(21) do not yet
   provide a bound on boundary source spectral density. The connected
   two-source comparison, finite-host estimates and retention inequality
   are also open. An argument purporting to close the task without these
   is a failed interpretation, not a result in this packet.
6. **Full geometry.** Check that no claim silently applies the fixed-record
   gap or old fixed-time buffer bound at times proportional to record size.

## Reproduction

From the repository root, run `python
research/reality_framework/grut_b2_scaling_analysis_20261018/verify.py` for
integrity. The frozen control command was:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python research/reality_framework/grut_b2_scaling_analysis_20261018/controls.py /tmp/b2_sa_review
```

It computes exact scalar fractions and three small full-matrix identity /
spectral-count controls. It computes NO time evolution or operating output
probability. Compare CONTROLS.json with the archived file; execution metrics
vary. Matrix arithmetic is ordinary double precision without outward
rounding; these controls do not certify the all-size proofs.

## Freeze chronology

* Parent B2-DX: `2d7207624d84e4aa24eac50abaebc112a7e0b267`.
* Charter and input pins: `8b1b9b8cd04d81b1331fd1d85ad5b35bdb0c873b`.
* Derivation, controls and source hashes: `1728dc630e476fae2534aa8856996e0cf7c609df`.
* Then one control execution and one byte-identical data replay.
* ROUTE_AND_DECISION.md (18) is a later analytical extension, explicitly
  separated from those frozen controls. No changed primary run or threshold.

No new primary operating experiment was launched. There is no L>=24
trajectory to reinterpret, no favorable-time search, and no B3.

## Literature boundary

Nachtergaele, Sims and Young, *Quasi-Locality Bounds for Quantum Lattice
Systems. Part I*, J. Math. Phys. 60, 061101 (2019),
https://arxiv.org/abs/1810.02428, is the primary context for propagation
upper bounds and quasi-locality. No theorem from it supplies our missing
outgoing probability, band-source weight, or formation guarantee. Our
finite-resource contour and Jacobi derivations are given explicitly for
this model. Recurrence literature and waveguide analogies are not imported
as diagnoses of the present family.

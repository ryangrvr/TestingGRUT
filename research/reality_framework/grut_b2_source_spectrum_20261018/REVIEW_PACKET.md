# Source proof review route

No independent reviewer has checked this packet. One assistant performed
all roles, using shared context and the inherited B1 physical catalog.

## Load-bearing checks

1. OBSERVABLES.md: fuel has scalar rank one ONLY in the one-resource source.
   Prove the operator G(t)I_2 and input-reference statement. Check the
   unnormalized comparator's joint event is a_L^2 and that the sufficient
   34% source condition maps to 11.458225% joint success, not 33.8%.
2. SOURCE_PROOF.md (7)-(9): rederive the two-component Rayleigh trial and
   the complement bound at k=pi. At the source root, check BOTH conversion
   directions and the counterterm: the residual must be d/2-J^2/d.
3. Equation (10): expand the full positive-supersolution identity, including
   all axes, fuel and charge configurations. Justify form-core closure and
   the no-eigenvector conclusion at threshold. It excludes no embedded or
   above-band point spectrum. Check the Weyl-sequence essential-band claim.
4. Equations (12)-(14): verify Cantelli/Chebyshev against the ORIGINAL fuel
   mean and variance. Recheck the bulk exact-graph extension used to place
   the central interval in the useful velocity range. This is an energy
   interval, not a source spectral-type or arrival theorem.
5. Equations (15)-(17): independently validate B1's original intertwiner
   residual and its half-line limit. Check the convex continuous-time upper
   bound, rather than sampled values. Check source fuel vs joint q=2, the
   already-frozen comparison, and the meaning of separated constituents.
6. SOURCE_VERDICT.md: the unresolved source spectral measure and observable
   kernel must remain unresolved. No finite pure-point measure can prove
   the proposed infinite-source total-variation comparison to a continuous
   reference by direct numerical agreement; their variation norm is 2.

## Reproduction and chronology

Parent `370098686ac5230dbf3a2619e793d11466775874`.
Charter/pins freeze `20a816e6af3c99f0baf302fda2d42cac5cd8b860`.
Proof/source freeze `f9b920870e58f55c713dbf7d0288f40a37fc0fff`.
Only then did the one primary control and one replay execute. Historical
proofs, comparison and law were not edited.

```sh
python research/reality_framework/grut_b2_source_spectrum_20261018/verify.py
```

Frozen source controls can be reproduced to a separate folder:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python research/reality_framework/grut_b2_source_spectrum_20261018/controls.py /tmp/b2_ss_review
```

There are no trajectories. K=8,12,16 matrices are complete finite source
compressions, with the artificial far-end counterterm removed. Every
relative position inside them is retained. R=32 is only a bulk arithmetic
control, not an inserted source cutoff. Floating positivity/eigenvalues
are not outward-rounded proof. Exact arithmetic applies to the rational
constants and interval endpoint comparisons only.

## Literature boundary

The user's probability-distance link was inspected:
https://qubit.guide/12.2-how-far-apart-are-two-quantum-states . Its elementary
norm-to-probability argument is reproduced directly in OBSERVABLES.md and
extended to the actual comparator with norm <=1. It supplies no source
dynamics or independence assumption. The positive-supersolution identity
is derived explicitly here; this is a model-specific application of known
spectral techniques, not a new physical principle or external validation.

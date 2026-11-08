# M0-01: the bounded finite-resolution matching gate missed

Frozen equations, implementation and analysis: commit `2b13287ae398d45bed7b4dcdcbb083747aaa4a00`, published before production. Executed from a clean detached checkout. [Protocol](PROTOCOL_M0_01.md), [machine card](protocol_m0_01.json), [raw analysis](raw_m0_01/analysis.json), [provenance](raw_m0_01/provenance.json).

The independent bosonized interval reference and the lattice solver diagnostics pass their registered checks. The finite lattice ladder does not reach the registered history window. This is a missed finite-resolution gate, not evidence of new physics or a failure of the continuum Schwinger model.

| Sites | Single-time TV | Times (0,1) TV | Times (0,2) TV | Times (0,1,2) TV |
|---:|---:|---:|---:|---:|
| 10 | 0.00315683 | 0.01119042 | 0.00970132 | 0.02010060 |
| 14 | 0.00238907 | 0.00808398 | 0.00739615 | 0.01425203 |
| 18 | 0.00191170 | 0.00622107 | 0.00598539 | 0.01110188 |
| 20 | 0.00173782 | 0.00558515 | 0.00545678 | 0.01001135 |

Every primary history improves monotonically along the fine three points, but every finest-point distance exceeds `0.001`. The three-time 18-to-20 step is `0.00109755`, also above its `0.00025` window. No extrapolated result is substituted for the frozen gate.

## What is controlled

- The independent direct-oscillator verification differs from the Gaussian word calculation by TV `8.41e-17`; the occupation-cutoff change is `6.94e-18`.
- The analytic three-time mode-tail bound is `5.50e-7`. Doubling continuum modes changes the evaluated probabilities by at most `7.76e-14`; tightening quadrature gives no resolved change at the reported precision. The analytic tail bound, quadrature estimates and measured changes remain different quantities.
- The tighter lattice solver and split propagation at N=14 change histories by TV `4.35e-14`. All registered eigen-residual and probability-normalization diagnostics pass. These diagnostics are not a general certified continuum or eigenstate-error theorem.
- Both primary and bare-zero controls completed within the registered resource caps.

## Where the discrepancy appears

The continuum mean is `0.1272298963`, variance `0.0621634882` and first gap `0.8444489696`. At N=20 the improved-mass lattice mean is `0.1114196229`, variance `0.0639852682` and gap `0.7977399530`. Mean and gap errors therefore remain substantial even though the history algorithm is well verified.

The bare-zero N=20 control has a smaller gap error but a larger maximum history distance, `0.02824370`. A successful spectrum match alone does not certify histories. Setting the continuum mean incorrectly to zero changes its three-time distribution by TV `0.08081654`; a zero-mean reference would compare the wrong preparation.

The observed trend is consistent with an unresolved leading cutoff/boundary correction. That is a diagnosis to test, not a proven allocation of the error. In particular, the geometric site span `(N-1)a` and the long-wavelength free open-chain length `(N+1)a` differ by an order-a endpoint correction. Any alternative dictionary must be derived and separately frozen; it cannot repair this verdict retrospectively.

## Scope of the result

This is the first direct comparison in this packet to independently calculated continuum detector histories, including their nonzero mean and three-time disturbance. It is massless only. No massive-history bound, uniform equation-(14) theorem, causal implementation of the ideal extended instrument, or reliability certification for an evaluator of unfamiliar theories follows.

The outcome is **finite-resolution matching unresolved at the registered accuracy**. Both failed windows and the passing implementation checks are retained. Zenodo remains paused.

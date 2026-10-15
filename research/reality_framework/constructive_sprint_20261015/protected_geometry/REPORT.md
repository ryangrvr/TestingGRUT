# Supplied geometry and matter response in the same evolving state

The frozen CP-12 Hamiltonian preserves a supplied flux-carrying square while supporting nonzero, gauge-invariant matter-density response on that same evolving preparation. This is a conditional mechanism; the graph and initial flux are explicitly supplied.

All 15 labelled four-edge graphs on four vertices are included. After unforced wait 0.5, a phase pulse exp(−iθn₀) is followed by lag 0.5. Receiver 2 is at initial-host distance two. No loop filter is used.

| Flux L | Host probability | Distance-two susceptibility | Truncation error | Response from escaped host |
|---:|---:|---:|---:|---:|
| 1 | 0.999998106615914 | -0.00015675148284 | 1.34e-13 | -2.04e-08 |
| 2 | 0.999999999999996 | -0.000206033144051 | 1.61e-13 | -4.16e-16 |

The response is overwhelmingly carried by the original-host component. Joint host–density responses are saved separately, so rare changed-graph outcomes are not being mistaken for transmission through the retained host.

Cutoffs |E|≤6 and |E|≤10 agree within the analytically bounded outgoing error and a reported floating-point tolerance. Every retained Hamiltonian column was independently compared with the unregulated integer transition implementation; every omitted transition was retained in the leakage operator. Bounds are ordinary floating-point evaluations, not formal interval-arithmetic certificates. Reported intervals add a separate numerical allowance of 5×10⁻¹²; this exceeds the observed cutoff discrepancies but is not a certified roundoff bound. The table separately lists the analytic electric-truncation estimate.

Two exact results support the mechanism. A supplied host with initial edge fields (±L,∓L) obeys

\[1-p_g(T)\le Q\left[\frac{(\sqrt2|t|T)^{2L}}{(2L)!}\right]^2.\]

Without any flux assumption, a fixed host also obeys $\|U(T)\psi-e^{-iP_gHP_gT}\psi\|\le [\Gamma\sqrt{Q(M-Q)}/M]T$. For sparse hosts this error is O(T/√N). These are supplied-host preservation results; the original neutral coherent input occupies many hosts and is not covered by the fixed-host root-mean-square estimate.

Exact signed, unregulated Taylor calculations give graph-loss onset T⁸/259200 for L=1 and (216017/12906142470144000000)T²⁴ for L=2. Higher-flux onset certificates and the complete Gauss target catalog are in EXACT_ONSET.json.

The test does not establish a macroscopic causal cone, autonomous geometric formation, preferred dimension, or gravity. The four-vertex host has diameter two. Initial electric energy and localized matter are physical preparation resources and are not matched to the original neutral ray.

Reproduce:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python exact_onset.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python report.py
```

## Neutral-origin recorded-history extension

A separate executed prototype begins from the original neutral ray, evolves unforced to time5, then applies a complete graph-and-circulation diagnostic record before the same gauge-invariant source probe. Charged full-circulation records have probability 0.055022107666. At lag0.5, their conditional host survival is 0.999999997486, circulation survival 0.893146518185, and opposite-on-recorded-ring density susceptibility 0.000131168601283 with interval [0.00013101254184763706, 0.00013132466074685398].

The complete external record costs system energy +0.0654739727471 and measured-system entropy 4.0824009083 bits. All outcomes are retained; the rare charged branch is analyzed conditionally. This joins formation and propagation in one finite recorded history while explicitly costing the diagnostic intervention. It is not an autonomous geometric phase.

See NEUTRAL_RECORD.md for exact symmetry handling, normalization, bounded energy-cost error control and numerical qualifications.

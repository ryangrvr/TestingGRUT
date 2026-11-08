# M1-01: interacting reference constructed; matching remains unresolved

Frozen source: `479eaed2f906e323dac491f8f308c0c05718d63b`. All calculations ran from a clean detached checkout, without refitting the mass, boundary conditions, detector or thresholds. [Protocol](PROTOCOL_M1_01.md), [analysis](raw_m1_01/analysis.json), [provenance](raw_m1_01/provenance.json), [events](raw_m1_01/events.json).

The independent bosonized calculation and lattice calculation give similar massive detector histories, but the registered reference-cutoff, overflow and final matching gates **miss**. The remaining reference error is large enough that this is not an independently certified interacting continuum bridge.

| Quantity | Measured | Frozen requirement | Verdict |
|---|---:|---:|---|
| Maximum N=20 versus E=12 history TV | 0.00247220 | <=0.001 | Miss |
| Bosonic E=10-to-12 history change | 0.00140056 | <=0.0001 | Miss |
| Maximum E=12 instrument overflow | 0.00126337 | <=0.00001 | Miss |
| Quadrature 128-to-256 history change at E=10 | 2.87e-14 | <=0.000001 | Pass |
| Lattice N=18-to-20 history change | 0.00011473 | <=0.00025 | Pass |
| Tightened lattice solver change | 1.32e-14 | <0.00000001 | Pass |

These numbers are absolute total-variation distances, not relative errors on individual rare outcomes. A small final step is a refinement diagnostic, not an absolute error bound.

## What was independently implemented

The reference uses a Dirichlet bosonic Fock basis and the normally ordered cosine Hamiltonian, without importing any lattice Hamiltonian, state, spectrum or covariance. Its normal-ordering factor agrees with a separate integral representation to 1.33e-15; the vertex matrix agrees with a large independently exponentiated oscillator to 2.50e-16. Its massless histories pass the frozen analytic-Gaussian comparison. All registered Hermiticity, eigen-residual, contraction and probability-normalization diagnostics pass. These checks establish numerical consistency within their scope, not a bound on the omitted interacting high-energy states.

The lattice comparator explicitly reuses the previously reviewed ED Hamiltonian. Its massive preparation is computed afresh with the declared bulk mass improvement and ghost-endpoint map. This is independent of the bosonic implementation, but not an independent implementation of the old lattice solver.

| Free-energy cutoff | Basis dimension | Massive finite-interval gap |
|---:|---:|---:|
| 4 | 21 | 1.63296256 |
| 6 | 72 | 1.62343276 |
| 8 | 207 | 1.61489681 |
| 10 | 534 | 1.61124091 |
| 12 | 1280 | 1.61018243 |

The gap improves with the cutoff, but that does not certify the histories. The absolute ground energies are cutoff-dependent and are not reported as renormalized vacuum energies. Normal ordering and subtracting a common free vacuum constant do not remove every Hamiltonian-truncation error.

## Why the instrument convention matters

The code uses the matrix elements of the uncut Weyl operator and completes its compression with an explicit overflow channel. It does not exponentiate a truncated field and call that the same measurement. The retained probabilities plus overflow sum to one. The 0.00126337 overflow is therefore visible rather than hidden by renormalizing the ordinary outcomes.

Nor is overflow alone the complete error: discarded Hamiltonian couplings and coherent returns are separate. The E=10-to-12 variation is itself above the declared tolerance. Increasing only arithmetic precision or spatial quadrature would not resolve that discrepancy.

## Physical conclusion

This pass has produced a fully explicit interacting reference and a concrete comparison of the same intended detector histories. It has not established their agreement to the registered accuracy. The two formulations have not been shown physically inconsistent: the reference truncation and finite-interval dictionary remain unresolved at that accuracy.

The bounded M1-01 run ends here. A reduced description of this same E=12 finite parent can be studied separately to identify which coherences its histories require, but it cannot certify the parent's continuum limit or replace this miss. No new GRUT law, experimental confirmation, universal equation-(14) statement or release follows. Zenodo remains paused.

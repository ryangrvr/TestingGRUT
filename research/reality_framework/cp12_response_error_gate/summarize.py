"""Build the scientific report from saved response records."""
from pathlib import Path
import json

P=Path(__file__).resolve().parent

def main():
    a=json.loads((P/'N6_B6_S64.json').read_text())
    b=json.loads((P/'N5_B10_S64.json').read_text())
    checks=json.loads((P/'CHECKS.json').read_text())
    rs=a['responses'];passed=all(r['sign_resolved'] for r in rs)
    table='\n'.join('| '+r['channel']+' | '+f"{r['projected_response']:+.10e}"+' | '+f"{r['observable_absolute_error_upper']:.6e}"+' | '+str(r['sign_resolved'])+' |' for r in rs)
    old=rs[1]['global_absolute_error_upper'];new=rs[1]['observable_absolute_error_upper']
    heading='SIX-EVENT DENSITY RESPONSES RESOLVED' if passed else 'SIX-EVENT RESPONSE ERROR GATE REMAINS PARTLY UNRESOLVED'
    report=f'''# CP-12 source response: sharper observable error control

Status: **{heading}; SPACETIME FORMATION NOT ESTABLISHED.**

Parent: `4570ed1c2ac2f8dab3f49ddb99c6224865eb1527`. This is an additive continuation of the unmasked formation benchmark. It does not change the frozen Hamiltonian, preparation, source, waiting/lag times, or electric budget. The new bound is our analytical result with independent numerical controls inside this packet, not an external mathematical audit or interval-arithmetic certificate.

## Result

For N=6, Q=12, the initial ray remains the neutral equal-amplitude superposition over every 12-edge graph. Mean energy is 6, energy variance 24, and initial min-entropy zero. It carries no topology mask, but its purity and specific boundary ray are supplied. The vertex source kick occurs at w=0.1 and the response is read after lag 0.1, in the frozen dimensionless CP-12 units.

| Detector | Projected susceptibility | Absolute error upper bound | Sign resolved |
|---|---:|---:|---|
{table}

These error bounds are derived analytically and evaluated in ordinary floating point. They are not statistical error bars or SI-unit predictions. The previous global state-error bound for the same distance-two channel was {old:.6e}; the new detector-specific bound is {new:.6e}, a factor of {old/new:.1f} smaller. The projected response itself is unchanged within numerical roundoff.

The five-event complete-graph control gives χ₁={b['responses'][0]['projected_response']:.10e}, with the new error upper bound {b['responses'][0]['observable_absolute_error_upper']:.6e}. Its distance-two detector is identically zero as an operator, since K₅ has no distance-two pairs. Our uniform global-norm bound deliberately still assigns a small conservative uncertainty to that zero detector; it need not be interpreted as an unknown physical channel.

## What changed in the calculation

The prior estimate bounded a measured bilinear response directly by the full wavefunction error. It included omitted amplitudes that a detector diagonal in electric-field space cannot directly read.

For a cutoff-commuting Hermitian detector O with global norm at most one, write e_x and e_y for the two forward state-error bounds, and d_x and d_y for backward residual bounds starting from detector-weighted projected output vectors. The new response-error upper bound is

    E_O = 2(e_y d_y + e_x d_x + e_x e_y).

DERIVATION.md proves this using Duhamel residuals and orthogonality of retained and omitted spaces. It does not assert that every observable enjoys a small error; backward errors must actually be computed. A hostile two-state example demonstrates why [O,P]=0 is essential.

Every outgoing full-law transition remains in the leakage matrix. The positive full rewiring diagonal is unchanged. The 147,153-vector signed source quotient and 180,850 outgoing orbits are the same as in the parent run. The numerical projection is still Σ E_l²≤6; the bound accounts for the unretained, unbounded-field dynamics.

The implementation batches independent backward probes and reuses exponential time grids. It uses the same parent midpoint Taylor upper sum. A dedicated control compares the grid and parent segment-by-segment implementations, including both endpoint vectors and upper sums. This improves computation reuse without asserting new physics. The largest three-column time grid here contains approximately 1.82 GB of complex entries before other solver memory; use adequate RAM.

## A matter-density response, with an important distance qualification

The first two detectors retain the parent's definition: density averaged over vertices lying at FINAL graph distance one or two from the source. Their response can include graph redistribution. Resolving them does not, by itself, isolate material transport through a geometry fixed at the intervention time.

We therefore also calculate the independent detector

    O_other = (1/(N−1)) Σ_(j≠0) n_j.

It contains no graph-distance projector. All six-event, 12-edge graphs are connected with diameter two, so O₁+O₂=O_other exactly on the full physical space. The code checks both the operator identity on the retained configurations and the corresponding response sum, and evaluates O_other's uncertainty directly with its own backward residuals.

Consequently, a resolved nonzero O_other response cannot be explained solely by relabelling the same receiver density among distance bins. It is a response of the mean matter population at vertices other than the source.

An exact attribution control also separates this result from the earlier cubic rewiring signal: with hopping removed, all remaining frozen terms commute with the source n₀. These detectors also commute with n₀, so every source-kick susceptibility is then identically zero for any initial state. Thus a nonzero mean-density response requires matter hopping. This operator control does not claim that hopping alone selects spatial topology or that every distance-bin effect is matter transport.

The source and all detectors are gauge invariant. The local phase kick is a diagnostic intervention applied to copies of the evolving state, not an autonomous preparation controller.

## What remains unestablished

The exact finite-size saturation obstruction from the parent still holds: in this sector every six-event graph has diameter two, B(1)=5, and B(r)=6 for r≥2. The geometric curves recorded with these source runs remain saturated. No three-dimensional volume-growth fit, photon identification, signal velocity, or causal cone is inferred.

The initial pure state remains a chosen boundary condition. Its abstract preparation accounting is documented in the parent; no physical state-selection principle has been added. Q/N=2 is also a chosen sector, not a dynamically selected connection density. No larger-size regulator-controlled spatial phase is demonstrated, and the inconclusive seven-event long-time pilot is not reclassified by these six-event results.

This result crosses a finite-time NUMERICAL RESPONSE precision gate. It does not cross the joint gate of autonomous preparation, persistent spatial organization, growing intrinsic dimension, and causal propagation. In particular, it cannot choose among all-to-all rewiring, available states, and missing state selection as the physical cause of a failure of that larger gate.

## Verification and reproduction

The verifier reruns the parent benchmark controls and the new tests. The new tests include {len(checks['dense'])} dense Hermitian examples with noncommuting dynamics, three complete small physical-matrix comparisons, the hostile noncommuting-projection example, and time-grid consistency. Actual direct-evolution response errors are checked against the derived upper bounds. The complete physical comparisons are independently built at half-leg cutoff two; they test the numerical implementation, while the unbounded-law error statement follows from the separate Duhamel proof.

Saved large-run records are checked for the bound arithmetic, the response-sum identity, global detector norm bounds, and numerical conservation. The verifier does not silently rerun those expensive records. The new proof is available for independent audit in DERIVATION.md, and the executable checks and declared settings are in README.md.
'''
    (P/'REPORT.md').write_text(report)

if __name__=='__main__':main()

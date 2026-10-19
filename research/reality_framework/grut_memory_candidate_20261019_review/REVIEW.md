# Review of GRUT-M1-W (commit 1997d20): numbers reproduce; the added principle does no work for retention

**19 October 2026. Branch `claude/m1-review-20261019`, on top of 1997d20. Reviewer: Claude (claude.ai), owner-directed, independent of the author. Additive; nothing in the reviewed packet is modified.**

## 1. Reproduction

`uniform_field_control.py` rebuilds H_M1 from the displayed formula by an independent tensor-product construction (not the packet's bit-table code) and re-runs the two core experiments at L = 3, 4, 5. Every reported core number is reproduced to the printed precision:

| quantity (L=4) | packet | this review |
|---|---:|---:|
| min retention success, t ≤ 10 | 0.987510 | 0.987510 |
| retention success at t = 40 | 0.935065 | 0.935065 |
| best matter readout | 0.778468 at t = 0.58 | 0.778468 at t = 0.58 |
| L=3 min retention / readout | 0.994999 / 0.806447 at 1.20 | 0.994999 / 0.806447 at 1.20 |
| L=5 min retention / readout | 0.989132 / 0.778460 | 0.989132 / 0.778460 |

The Hamiltonian is Hermitian as displayed (F_i acts on the neighbours of i and commutes with τ_i^x). The writer is not re-run here; its mechanism is checked analytically in §3.

## 2. The missing control, and what it shows

The packet's "unfacilitated" control sets F_i ≡ 1 (transverse field Δ = 1 on every record spin). That is not the right baseline for the claim that *kinetic saturation protects records*: a transverse field of order the Ising coupling destroys any ferromagnetic memory, facilitated or not. The right baseline keeps the same floor amplitude everywhere and removes only the facilitation: **F_i ≡ ε**, i.e. a uniform transverse field εΔ = 0.05 with no neighbour dependence. Everything else (matter chain, α-feedback, K, parameters, states, observables, time grids) is identical.

Result (`control_results.txt`):

| L | model | min success t ≤ 10 | success at t = 40 | min over [0, 40] | best readout |
|---|---|---:|---:|---:|---:|
| 3 | M1 (facilitated) | 0.994999 | 0.985006 | 0.984591 | 0.806447 |
| 3 | uniform field 0.05 | **0.998738** | **0.994692** | **0.994266** | 0.806626 |
| 4 | M1 (facilitated) | 0.987510 | 0.935065 | 0.935065 | 0.778468 |
| 4 | uniform field 0.05 | **0.999055** | **0.998458** | **0.998397** | 0.778483 |
| 5 | M1 (facilitated) | 0.989132 | 0.940399 | 0.940399 | 0.778460 |
| 5 | uniform field 0.05 | **0.998479** | **0.997517** | **0.995717** | 0.778475 |

The control retains the record **better** than the candidate at every size and every time, by a factor of 3 (L = 3) to 40 (L = 4) in the error at t = 40 (L = 4: 0.15 % vs 6.5 %). The matter readout is unchanged to four decimals, because readout depends only on the α τ^z σ^z σ^z feedback, which both models share.

**Reading.** In the demonstrated regime the memory is the ferromagnetic gap of −K τ^z τ^z against a small transverse field; that is a zero-temperature transverse-field Ising repetition code, and it works for the same reason in both models. The facilitation term F_i → 1 at a mismatched neighbourhood is inactive on an ordered record (every neighbourhood is ordered), so it cannot help retention at first order; and once a single record spin leaks, it turns both neighbours' flip amplitude from 0.05 to 1, which converts an off-resonant single-flip leakage into resonant spreading (an East-model cascade). That is why the candidate loses more than the control. Kinetic saturation, as implemented, is a **liability** for retention of a prepared record, not its mechanism. The packet's own "unfacilitated" comparison (0.62) is a comparison with a model that has no memory for an unrelated reason.

The packet already notes that facilitation is also a barrier to formation (domain-wall conservation at ε = 0). With this control, the term is seen to be neutral-to-harmful on both sides of the write/protect trade-off it was introduced to manage. What the construction demonstrates is therefore: *a gapped classical repetition code coupled to a matter chain through a τ-dependent bond, read out through the matter.* That is a correct and clean finite model; it is not evidence for the "state-dependent kinetic saturation" principle.

## 3. The writer

Parameters: v = p₀/M = 200, gate half-width w = 0.3, so the carrier crosses a gate in 2w/v = 0.003, three orders of magnitude below every core time scale (J = K = Δ = 1). The pulse area is Ω ∫ f_i dt = Ω w / v = π/2 by construction, so the conditional coupling applies e^{−i(π/2) τ_i^x} = −i τ_i^x: an exact controlled flip, with the register frozen during passage. The ≥ 0.99968 per-cell fidelity is the residual from finite packet width and the small core evolution during 0.003. The energy transfer (core 1.55 → −0.44, carrier +1.99) is the register's energy change under the flip, taken from the carrier's kinetic energy, as it must be in a closed unitary model.

So the writer is an impulsive flying-qubit controlled-NOT on each cell. It is a legitimate closed-system apparatus, and the packet says correctly that it is "a programmed writer with prepared input and fuel". It does not test anything about the memory law: the same particle writes the uniform-field control equally well.

## 4. Where this leaves M1 relative to the rest of the record

- The companion packet `toric4d_memory_vs_scrambling_20261019` (branch `claude/toric4d-memory-vs-scrambling-20261019`) tested the same qualitative statement, "records persist while local observables spread", in the one model where the memory is *proven* (4D toric code), at finite temperature, with size scaling, and found it true, and found that it gives no non-scrambling. M1 at L = 3–5 and T = 0 is consistent with that but adds nothing to it: a 1D Ising ring has no finite-temperature memory (Peierls), so the regime in which M1 works is the regime in which any gapped ferromagnet works.
- The one new ingredient M1 proposed, kinetic saturation, is shown here to be inert or harmful for the record it was meant to protect. The honest classification of M1-W is: **a correct finite model of a repetition code with matter feedback and a flying-qubit writer; the claimed principle is not supported by its own numbers once the right control is run.**
- The packet's own boundary statements (no gravity, no Lyapunov exponent, no priority, no robustness) stand and are accurate.

## 5. Suggested amendment to the M1 packet

Add the F_i ≡ ε control to `calculate.py` and `summary.json`, and change the headline from "kinetic saturation protects records" to "a gapped record couples to matter through the bond; the facilitation term is not responsible for retention in the demonstrated regime". If the author wants a regime where facilitation helps, it has to be one where mismatched neighbourhoods are *supposed* to relax fast while ordered ones hold — i.e. error correction by cascade — and that has to be shown against the uniform-field baseline, not against Δ = 1.

Reproduction: `python3 uniform_field_control.py` (about one minute; numpy, scipy).

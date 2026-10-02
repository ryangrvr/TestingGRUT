# ARROW ORIGIN LEDGER

Every claimed effective arrow (subsystem, record, coarse) gets a row. "Typical" must name its measure.

| ID | Hamiltonian / law | TPS Σ | global initial state | environment marginal | initial S–E correlation | bath freshness | measure / typicality rule | coarse-graining A | time orientation | subsystem arrow? | record arrow? | exact info preserved? | price |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AR-01 | any finite closed unitary H | any | **every** admissible state | — | — | — | none | any continuous A | — | **NO** (Theorem 1 / Cor. 1) | NO | yes | — (no-go) |
| AR-02 | mixed-field Ising, n = 10 | computational, S = 2 qubits | real product (tilt 1.4) | pure product | 0 | — | none | reduced | forward and backward symmetric | yes, finite window, two-sided (Janus) | — | yes | H_CORR-PRICED, Σ-PRICED, A-PRICED, TIME-ORIENTATION-PRICED |
| AR-03 | same | same | time-reversed Kψ(t*) | entangled | high | — | none | reduced | — | **anti-arrow** (1.238 → 0) | — | yes | shows the every-state no-go |
| AR-04 | same | same | Haar / energy-shell typical | ~maximal | ~maximal | — | **Haar (shell)** | reduced | none (P(inc) = 0.499 ± 0.013) | **no** (equilibrium only) | no | yes | MEASURE-PRICED (equilibrium ≠ arrow) |
| AR-05 | mixed-field Ising, n = 8 | S = qubit 0 | product \|0⟩⊗𝟙/128 | **maximally mixed** | 0 | — | none | reduced | forward | yes | — | yes | H_CORR-PRICED (independence), **not** low-entropy |
| AR-06 | same | same | correlated Kρ(τ)K | 6.902 bits | I = 0.857 | — | none | reduced | — | **anti-arrow** (0.955 → 0) | — | yes | H_CORR (correlations reverse the arrow) |
| AR-07 | same | same | product of AR-06 marginals | 6.902 bits | 0 | — | none | reduced | forward | yes (relaxes) | — | yes | H_CORR-PRICED |
| AR-08 | partial-swap collisions | S \| ancillas | S ⊗ fresh ancillas | pure / thermal / **max-mixed** | 0 | **fresh** | none | system marginal | forward | yes, Markovian | **yes** (0.225 – 0.324) | yes (global) | FRESHNESS-PRICED, H_CORR-PRICED |
| AR-09 | same | same | S ⊗ correlated ancillas (classical / GHZ) | max-mixed marginals | bath-internal correlations | fresh | none | system marginal | forward | entropy yes; **non-Markovian** | yes | yes | H_CORR (bath correlations break the semigroup) |
| AR-10 | same | same | S ⊗ reused bath (M = 1, 2, 4, 7) | max-mixed | builds up | **not fresh** | none | system marginal | none | **no** (recoherence; M = 1 repurifies) | records erased | yes | FRESHNESS-PRICED |
| AR-11 | Σ g_k Z_S Y_k | S \| 6 env | \|+⟩\|0…⟩ | pure product | 0 | — | none | fragments | forward | — | **yes** (redundancy 0 → 6) | yes | H_CORR, Σ, A (fragments) |
| AR-12 | same | same | Θψ(t_f) | — | high | — | none | fragments | — | — | **records un-form** (6 → 0) | yes | RECORD ARROW H-PRICED |
| AR-13 | mixed-field Ising, n = 10 | **3 frames** | the same ψ_prod(4) | frame-dependent | frame-dependent (0 in the Householder frame) | — | none | reduced | frame-dependent | frame-dependent | — | yes | **Σ-PRICED** (low correlation not TPS-invariant) |
| AR-14 | same | computational | mean-field product state (fixed by (H, Σ) mod reflection) | product | 0 | — | none | reduced | **two-sided** | Janus | — | yes | STRUCTURALLY SELECTED given (H, Σ, selector functional); TIME-ORIENTATION NOT SELECTED |

# Does a passive self-correcting memory forbid local scrambling? (pre-registered test)

**19 October 2026. Branch `claude/toric4d-memory-vs-scrambling-20261019`, additive. Author: Claude (claude.ai), owner-directed. Status: PRE-REGISTERED, NOT YET RUN at the time of this commit.** Frozen CP-12, Candidate 1, GRUT III and the Stage-3 charter are unchanged.

## 1. What is being tested, and why

The sharpest surviving reading of GRUT's founding intuition (closed whole; memory; saturation into rigid structure; geometry is the structure) is:

> the late-universe state is a passive self-correcting quantum memory whose rigidity is geometry, and therefore the fundamental dynamics does not scramble (λ_L < 2πT).

The literature brief of 19 Oct 2026 (research report, this conversation) found that the second clause is contradicted by every published computation at the de Sitter horizon (Aalsma–Shiu 2020 λ_L = 2π/β; Kolchmeyer–Liu 2024 and Narovlansky 2025 4π/β), and that the first and second clauses are joined only by an unstated **bridge claim**:

> **B:** a passive self-correcting memory is necessarily non-scrambling / non-thermalizing in its local degrees of freedom.

AdS/CFT is a standing counterexample class (a code that is maximally chaotic), but no direct test of B in a proven passive memory was found. This packet tests B in the only translation-invariant model that is *proven* to be a passive self-correcting memory: the 4D toric code (Dennis–Kitaev–Landahl–Preskill 2002; Alicki–Horodecki–Horodecki–Horodecki 2010). If local information thermalizes and spreads at the same temperature and size where the logical bit survives, then "memory" and "non-scrambling" are independent properties, B is false in the one model where it can be checked, and the unification that motivated the sharp reading is gone.

This is a classical stochastic test. Its scope and its limits are stated in §5 before running.

## 2. Models

**4D toric code, Z-error sector.** Qubits on the 2-cells (faces) of an L⁴ hypercubic torus, 6L⁴ faces. Edge stabilizers A_e = ∏_{f∋e} X_f (6 faces per edge). A Z-error configuration is a 2-chain e of faces; it violates A_e on the edges of ∂e. Hamiltonian restricted to this sector: E = 2·|∂e| (each violated stabilizer costs 2). A logical error is a membrane spanning two torus directions. Readout: logical X̄ on plane (01) is the set M_X = {(01)-faces at x₀ = x₁ = 0}; the stored bit is flipped iff |e ∩ M_X| is odd after decoding. The thermal model is the 4D Z₂ **2-form** gauge theory (variables on faces, constraints on their boundary edges), which is dual to the **4D Ising model**; the duality tanh β* = e^{−2K} with K_c(4D Ising) ≈ 0.1497 gives β_c ≈ 0.954. *(Pre-registration amendment, before any memory run: the first draft of this file said β_c ≈ 0.44, the 1-form gauge theory value; a smoke-test energy scan, kept as `results/energy_scan_prelim.json`, showed the system still disordered at β = 0.6 and a transition between β = 0.85 and 0.95 at L = 3 and 4, consistent with the 4D-Ising duality. The β grid below was changed accordingly; the readings in §4 were not.)*

**2D toric code control.** Qubits on edges of an L² torus. Z-errors are strings; violated vertex stabilizers are their endpoints; E = 2·#endpoints. No ordered phase (T_c = 0); the memory time is expected L-independent at fixed T. Readout by minimum-weight perfect matching (pymatching).

**Dynamics (both).** Single-site Metropolis at inverse temperature β: pick a random face/edge, flip with probability min(1, e^{−βΔE}). Local, Markovian, detailed-balance; a Glauber/Davies-type bath. Time unit = one sweep (one attempted flip per site).

**4D decoder.** Zero-temperature local cooling: Metropolis with ΔE<0 accepted and ΔE=0 accepted with probability ½, for L² sweeps. Contractible membranes have boundary tension and erode (the single-plane restriction is the 2D Ising model, whose T=0 droplets shrink); a spanning sheet has no boundary and survives. This is a local decoder, not minimum-weight; its failure mode (a large contractible membrane wandering into a wrapping one) is suppressed at low T and checked by the β>β_c vs β<β_c contrast.

**Quantum illustration (not decisive).** Exact 2D toric code on a 3×3 torus (18 qubits) with H = −ΣA_s − ΣB_p − h Σ(X_i+Z_i), h = 0.3. Infinite-temperature OTOC C(t) = 1 − Re⟨W(t)VW(t)V⟩ for local Paulis, stochastic trace. This shows only that a perturbed stabilizer code scrambles; 2D is not self-correcting, so it does not test B.

## 3. Observables

| ID | Observable | Meaning |
|---|---|---|
| O1 | P_log(t) = ⟨decoded logical bit(t) · bit(0)⟩, from the zero-error state; τ_log from P_log = e^{−t/τ_log} | memory time |
| O2 | Equilibrium autocorrelation of a single face spin and of a single edge syndrome, C_loc(t); τ_loc | local thermalization time |
| O3 | Damage spreading: two equilibrium copies differing by one flipped face, identical random stream; Hamming distance D(t) | local information spreading (update-rule dependent; secondary) |
| O4 | ⟨E⟩/N_faces versus β, L=4 | locate the transition |
| O5 | 2D control: τ_log(L, β) by MWPM readout | method control |
| Q1 | OTOC C(t) in the 18-qubit perturbed 2D toric code; local vs logical infinite-T autocorrelation | illustration only |

Grid: 4D L ∈ {3,4,5}, β ∈ {0.60, 0.90, 1.00, 1.10} (one well below β_c ≈ 0.95, one at the rounded transition, two above); 2D L ∈ {4,6,8,12}, same β; 16 trajectories per point; readout on a logarithmic time grid to t_max = 10⁵ sweeps; where P_log has not decayed by t_max a lower bound on τ_log is reported, not a fit. At β ≥ 1.0 the 4D memory time at L = 4, 5 is expected to exceed t_max; R1 then reads "τ_log(L=3) finite and τ_log(L≥4) > t_max > τ_log(L=3)" as an increase.

## 4. Pre-registered readings (written before running)

**R1 — method control.** At β > β_c, τ_log(4D) increases with L (expected roughly exponentially); at β < β_c it is O(1) for all L. τ_log(2D) is L-independent at every β. If τ_log(4D) does *not* increase with L at reachable sizes, the packet is INCONCLUSIVE on finite-size grounds and stops; nothing is inferred about B.

**R2 — the test of B.** At the same (β, L) where R1 shows memory, O2 gives τ_loc = O(1), independent of L, and O3 gives D(t) that does not stay at its initial value (it either grows toward an extensive value or heals to zero; either way the local record is not retained as a local record).
- Outcome R2-a (expected): τ_loc = O(1) and L-independent while τ_log grows with L → **local degrees of freedom thermalize while the logical bit persists; B is FALSE in the 4D toric code.** The two clauses of the sharp reading are independent. The non-scrambling clause stands alone and is already contradicted at the dS horizon by the literature. Verdict: the sharp reading loses its unifying content; retire it.
- Outcome R2-b (would be new): τ_loc grows with L, or D(t) stays bounded away from both 0 and the extensive value, at β > β_c only → the memory phase suppresses local thermalization in this model → B is ALIVE in a classical proxy → escalate to a quantum OTOC computation in a large-N self-correcting model before any further claim. No GRUT claim is made on R2-b alone.

**R3 — illustration.** Q1 shows OTOC growth and saturation under perturbation in a stabilizer code. Any outcome of Q1 is non-decisive for B.

**Reading rule.** The verdict is read off R1 and R2 by the table above. No reading is changed after the numbers are seen. If a script fails or a fit is ambiguous, that is reported as such.

## 5. Scope and what is not tested

- This is a classical stochastic (Davies-type) model of the thermal dynamics of a stabilizer code. It tests local thermalization and classical information spreading in a passive memory. It does not compute a quantum OTOC or a Lyapunov exponent in the 4D code (infeasible here: 96 qubits at L=2).
- The open question "can a passive self-correcting phase carry a Lorentz-invariant spin-2 low-energy sector?" is not touched.
- Nothing here concerns de Sitter space; the dS chaos results are taken from the literature brief.
- Damage spreading (O3) depends on the update rule and is reported as secondary evidence only; O2 (autocorrelation) carries R2.
- Priority: thermal stability of the 4D toric code is DKLP 2002 / Alicki et al. 2010; the 2D memory-time result is Alicki–Fannes–Horodecki; damage spreading as a chaos proxy is Derrida–Weisbuch 1987 (classical). The only thing this packet adds is the *joint* measurement of memory and local thermalization in the same trajectories, as a direct test of B.

## 6. Reproduction

`python code/run_all.py` writes `results/*.json` and `results/SUMMARY.md`. Requirements in `code/requirements.txt` (numpy, scipy, numba, pymatching). Seeds are fixed in the scripts.

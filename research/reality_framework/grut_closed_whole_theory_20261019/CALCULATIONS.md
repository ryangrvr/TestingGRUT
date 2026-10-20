# Calculations: one closed whole carried from preparation to measurement, and the compatibility check

19 October 2026. Packet `grut_closed_whole_theory_20261019`. Everything in this file is
standard quantum mechanics evaluated exactly in a **supplied** finite model. The purpose is
to make the working theory's claims about closure calculable and to show, with numbers,
what closure forces and what it leaves free. No fitted inputs, no new physics.

Code: `code/closed_whole_demo.py` (calculation 1) and `code/marginal_compatibility.py`
(calculation 2). Outputs: `results/closed_whole_demo.json`, `results/marginal_compatibility.json`.
Reproduce with `python3 code/closed_whole_demo.py` (about 1 minute, numpy + scipy) and
`python3 code/marginal_compatibility.py` (under a second, cvxpy + CLARABEL). Script hashes and
library versions are recorded in each JSON under `provenance`. ħ = 1 throughout.

## Calculation 1: the closed whole Ω = A ∪ Ā

### Theory datum (everything that had to be supplied)

| Item | Supplied value |
|---|---|
| Hilbert space | 8 qubits, dim 256 |
| Access structure | A = qubit 0 (the "system"); Ā = qubits 1–7 (its bath = the rest) |
| H_A | (ω₀/2) Z₀, ω₀ = 1 |
| H_Ā | XXZ chain on qubits 1–7: Σ_k [J(X_kX_{k+1} + Y_kY_{k+1}) + J_z Z_kZ_{k+1}] + Σ_k h_k Z_k, J = 1, J_z = 0.5, h = (0.3, −0.7, 0.5, −0.2, 0.8, −0.4, 0.6) |
| H_{AĀ} | g X₀X₁, g = 0.5 |
| Preparation | ρ_Ω(0) = |+⟩⟨+|_A ⊗ e^{−βH_Ā}/Z, β = 1 (a mixed global state; closure does not require purity) |
| Instruments | projective Z and X measurements on A; Born rule p = Tr(ρ_A E) |
| External time | the parameter t of U(t) = e^{−iHt}; grid t ∈ [0, 40], step 0.1 |

The complement is a disordered XXZ chain with no known integrability; H_Ā conserves the total
bath Z and the full H conserves the Z-parity of all eight qubits. Nothing below needs
non-integrability. The global law is ρ_Ω(t) = U(t)ρ_Ω(0)U(t)†; the local state is ρ_A(t) = Tr_Ā ρ_Ω(t).

### S1. Preparation → evolution → measurement, with accounting

| t | E_tot | E_A | E_Ā | E_int | S(ρ_A) | I(A:Ā) | p(Z=+1) | p(X=+1) |
|---|---|---|---|---|---|---|---|---|
| 0 | −10.228789 | 0.000000 | −10.228789 | 0.000000 | 0.000000 | 0.000000 | 0.500000 | 1.000000 |
| 1 | −10.228789 | −0.037419 | −10.018595 | −0.172774 | 0.135553 | 0.185858 | 0.462581 | 0.784788 |
| 2 | −10.228789 | −0.126915 | −9.933866 | −0.168008 | 0.335708 | 0.445762 | 0.373085 | 0.392354 |
| 5 | −10.228789 | −0.126617 | −9.928410 | −0.173762 | 0.355440 | 0.462109 | 0.373383 | 0.337724 |
| 10 | −10.228789 | −0.143593 | −9.882558 | −0.202637 | 0.392871 | 0.453948 | 0.356407 | 0.212817 |
| 20 | −10.228789 | −0.056339 | −10.011149 | −0.161301 | 0.333707 | 0.423519 | 0.443661 | 0.598989 |
| 30 | −10.228789 | −0.121234 | −9.917211 | −0.190343 | 0.448388 | 0.553655 | 0.378766 | 0.619657 |
| 40 | −10.228789 | −0.098988 | −9.983861 | −0.145940 | 0.419838 | 0.483332 | 0.401012 | 0.178277 |

Invariants over the whole grid: max |ΔE_tot| = 7.1·10⁻¹⁵; max |ΔS(ρ_Ω)| = 4.0·10⁻¹⁵ (S(ρ_Ω) =
1.260931); max |Δ Tr ρ_Ω²| = 6.7·10⁻¹⁶ (Tr ρ_Ω² = 0.489575). Local Born probabilities computed
from ρ_A agree with the global expression Tr(ρ_Ω (E⊗1)) to 2·10⁻¹⁶. The local energy E_A and
the interaction energy E_int are **not** conserved separately; only their sum with E_Ā is.
S(ρ_A) stays below ln 2 = 0.693 and first exceeds 0.4 at t = 8.8.

Exact reduced equation, checked at t = 3 by central differences (dt = 10⁻⁴):

    dρ_A/dt = −i Tr_Ā [H_A + H_{AĀ}, ρ_Ω]      (|lhs − rhs| = 1.9·10⁻¹⁰),

and the H_Ā commutator traces to zero identically (3·10⁻¹⁶). The right-hand side is a
functional of the **global** state, not of ρ_A: that is the whole content of "bath = rest".

### S2. Closure does not force monotone relaxation

S(ρ_A) is not monotone: it has 13 local maxima and 26 sign changes of dS/dt on the grid. The
first local maximum, 0.3383 at t = 2.1, is followed by a drop to 0.2353 at t = 3.6; the global
maximum 0.4893 occurs at t = 37.5 and S falls to 0.4042 by t = 39.5. Information flows back
from the bath into A: the trace distance D(t) between the reduced states of two preparations
that differ only on A (|+⟩ vs |−⟩, same bath state) starts at 1, reaches its minimum 0.5418 at
t = 37.6, and accumulates **1.2035** in positive increments on the 0.1 grid (a lower bound on
the continuous Breuer–Laine–Piilo measure; largest single step 0.0145). A CP-divisible reduced
dynamics would make D non-increasing, so the backflow certifies memory. The converse does not
hold: CP-divisible dynamics can still carry memory, and Markovianity is a property of the
multi-time process (Pollock et al., Phys. Rev. Lett. 120, 040405 (2018), primary text inspected).

### S3. The reduced dynamics is derived, and it is state-dependent

For product preparations with a fixed bath state, ρ_A(0) ↦ ρ_A(t) is a linear CPTP map Φ_t
(min Choi eigenvalue over sampled t: 3·10⁻¹¹; trace preservation to 3·10⁻¹⁵). Its exact
time-convolutionless generator K(t) = Φ̇_t Φ_t⁻¹, which exists wherever Φ_t is invertible
(4×4 in the Pauli basis; condition numbers 1.1–2.3 at the six sample times, at most 2.86 over
the full 0.1 grid), is time dependent and depends on which bath state was supplied:

| t | ‖K(t) − K(0.5)‖, thermal bath | ‖K_thermal − K_ground‖ | ‖K_thermal − K_{T=∞}‖ |
|---|---|---|---|
| 1 | 0.1598 | 0.0758 | 0.4201 |
| 2 | 0.5202 | 0.1001 | 0.6511 |
| 5 | 0.9689 | 0.1084 | 0.9237 |
| 10 | 0.6789 | 0.0688 | 0.4313 |
| 20 | 0.6010 | 0.0651 | 0.2669 |

There is no fixed Lindbladian, and the kernel is downstream of the supplied bath state. The
intermediate maps V(t+0.5, t) = Φ_{t+0.5}Φ_t⁻¹ have negative Choi eigenvalues at every
sampled t (−0.0038 at t = 0.5 down to **−0.1590** at t = 5): the derived reduced dynamics is
not CP-divisible.

### S4. Identical marginals, different futures

Two global preparations ρ₁ = ρ_A ⊗ ρ_Ā and ρ₂ = ρ₁ + εχ with Tr_A χ = 0 = Tr_Ā χ have the
**same** ρ_A(0) and the **same** ρ_Ā(0) (verified to machine precision), ρ₂ ≥ 0, and yet:

| Case | ε (= ½ ε_max) | ‖εχ‖₁ | max_t D(ρ_A^{(1)}(t), ρ_A^{(2)}(t)) | at t | A-factor of 𝓘(0), ‖−i Tr_Ā[H, εχ]‖ |
|---|---|---|---|---|---|
| T = ∞ bath, χ = Y₀X₁ | 1/512 | 0.500 | **0.1025** | 0.7 | 0.3536 |
| T = ∞ bath, χ = Z₀Z₁ (classical correlation) | 1/512 | 0.500 | 0.0097 | 14.1 | 0 |
| β = 1 bath, χ = Y₀ ⊗ √ρ_Ā (X₁ − ⟨X₁⟩) √ρ_Ā | 0.25 | 0.293 | 0.0544 | 0.7 | 0.1925 |
| T = ∞ bath, χ = Z^{⊗8} (conserved Z-parity, [H, χ] = 0) | 1/512 | 0.500 | **4·10⁻¹⁶** | n/a | 0 |

The last column is the A-factor of the Nakajima–Zwanzig inhomogeneity at t = 0 (𝓘(0) is this
operator tensored with ρ_ref). For χ = Y₀X₁ it is exact by hand:
[gX₀X₁, Y₀X₁] = 2ig Z₀ ⊗ 1_Ā, so −i Tr_Ā[H, εχ] = 2gε·2⁷ Z₀ = ¼ Z₀, Frobenius norm √2/4 =
0.35355, as computed. No map acting on ρ_A(0), nor on the pair (ρ_A(0), ρ_Ā(0)), reproduces
both futures: the reduced state alone does not specify an arbitrary preparation's future.
For the classical correlation Z₀Z₁ the first two orders of 𝓘 vanish at t = 0 and the difference
appears at third order in t, which is why it stays small. The fourth row is the converse
counterexample: the preparation is correlated (𝒬ρ(0) ≠ 0), yet [H, χ] = 0 makes 𝓘 vanish for
all t and the two futures coincide. The inhomogeneous term vanishes for preparations
ρ_A ⊗ ρ_ref, and not only for them. (Pechukas 1994; Štelmachovič–Bužek
2001; standard.)

### S5. Closure does not force a unique local law

Same H_A, same ρ_A(0) = |+⟩, same coupling form gX₀X₁, three different closed complements:
(a) the disordered XXZ chain above; (b) a single bath spin, H_Ā = (ω₀/2)Z₁; (c) a clean XX
chain (free fermions). Trace distances between the resulting ρ_A(t):

| t | a vs b | a vs c | b vs c | S(ρ_A): a / b / c |
|---|---|---|---|---|
| 1 | 0.0289 | 0.0126 | 0.0167 | 0.136 / 0.180 / 0.160 |
| 5 | 0.6046 | 0.2765 | 0.7107 | 0.355 / 0.168 / 0.441 |
| 10 | 0.6051 | 0.3527 | 0.4860 | 0.393 / 0.478 / 0.504 |
| 40 | 0.6541 | 0.7930 | 0.2328 | 0.420 / 0.410 / 0.149 |

With the single-spin complement and g′ = 0.75 all Bohr frequencies are multiples of 0.5, so the
closed whole returns exactly after T = 4π = 12.566: trace distance of ρ_Ω(T) to ρ_Ω(0) is
7.6·10⁻¹⁵, S(ρ_A) goes 0 → 0.4487 (T/2) → 0 (T). Relaxation is a property of the supplied
complement, not of closure.

### S6. Records and interventions: an externally timed write

Add a register qubit R (|0⟩), H′ = H ⊗ 1_R. At t₁ = 3 apply CNOT(A→R); evolve to t₂ = 6. The
CNOT is not generated by H′ (H′ conserves Z_R, so by THM 2.15 it can never write R): it is an
externally timed intervention, a supplied control and clock that open the whole at t₁ and
change its energy by 0.2076. This is the case the process-tensor framework describes (DEF 1.5
(a)). A strictly internal write needs a clock and controller inside Ω (DEF 1.5 (b)), as in the
program's M1 flying-qubit writer.

- Born probability p(Z_A = +1, t₁) = 0.354910; population of |0⟩_R after the write = 0.354910.
  This equality is an identity (a CNOT on a fresh |0⟩ is a Stinespring dilation of Z-dephasing
  on A), recorded as a consistency check: the record **is** the Born probability, stored in a
  subsystem.
- Joint readout at t₂: p(+,+) = 0.2155, p(+,−) = 0.1230, p(−,+) = 0.1394, p(−,−) = 0.5220.
- p(Z_A = +1, t₂): **0.338531** exact, with the record written; **0.449147** with no
  intervention; **0.288331** from composing one-time maps with the instrument
  (Φ_{t₂−t₁} ∘ dephase ∘ Φ_{t₁}, bath re-prepared at t₁); 0.338531 if A is dephased at t₁ with
  nothing kept (identical to the unitary record because the idle register commutes with the
  partial trace, again an identity). Multi-time statistics under intervention are properties of
  the process, not of Φ_t.
- Record persistence. With [H′, Z_R] = 0 the register is unchanged at all later times
  (7·10⁻¹⁵ at t = 10). Adding κX_R (κ = 0.3) rotates the record out of the fixed Z_R readout
  basis (trace distance 0.0858, 0.2273, 0.2505 at t = 4, 6, 10) but loses nothing: the optimal
  readout of the recorded bit from R stays at success 1.000000. Coupling the register to the
  bath, κX_RX₇, makes the bit genuinely less retrievable: optimal readout success 0.9607 (t = 4),
  0.9672 (t = 6), 0.9613 (t = 8), 0.9700 (t = 10). Whether a record persists is a property of the
  supplied dynamics, not a consequence of closure.

### S7. No arrow is forced

Forward evolution to T = 40 followed by U(−T) returns ρ_Ω(0) to 5·10⁻¹⁵. H is real in the
computational basis, so complex conjugation is a time reversal and ρ_Ω(T)* is an admissible
preparation of the same closed whole. Its forward evolution is the original run read backwards,
ρ_rev(t) = ρ_Ω(T−t)* exactly (deviation 5·10⁻¹⁵), so S(ρ_A) returns **non-monotonically** to
0 at t = 40 (sampled 0.4198 → 0.4005 → 0.4484 → 0.3337 → 0.3929 → 0.0000; it rises on 201 of 400
grid steps and reaches 0.4893 at t = 2.5). The same closed law carries a legal preparation from
local entropy 0.42 to 0; a low-entropy preparation, not closure, supplies the arrow.

### S8. What closure does force (checked for a pure global preparation)

With ρ_Ā(0) a pure eigenstate of H_Ā: S(ρ_A) = S(ρ_Ā) at all sampled times (deviation
3·10⁻¹⁶) and global purity constant to 10⁻¹⁴; the closest return to the initial state for t ≥ 1
within the grid is trace distance 0.401 at t = 11.1 (quasi-recurrence is a theorem in finite
dimension; its time is not short and no recurrence time is claimed). Subadditivity and
Araki–Lieb are trivial for a pure global state, so they are checked on the mixed S1 trajectory
(S(ρ_Ω) = 1.2609): the mutual information is ≥ 0 at all 401 times (minimum 1.8·10⁻¹⁵ at t = 0,
the product state) and the Araki–Lieb slack S_Ω − |S_A − S_Ā| is 0 at t = 0, 0.2115 at t = 3 and
0.2439 at t = 20, with no violation. All of these are theorems of standard quantum mechanics.

### S9. A unitarily equivariant state rule built from H alone is stationary

If a rule picks the global state from H alone and is equivariant under unitary conjugation,
F(UHU†) = UF(H)U†, then taking U = e^{iHs} gives [F(H), H] = 0: the selected state is
stationary, so its one-time marginals do not evolve when nothing intervenes. Checked for the
Gibbs state at β = 1 and the ground projector of the same eight-qubit H: ‖[F, H]‖ = 8.9·10⁻¹⁵ and
1.6·10⁻¹⁴, and ρ_A(t) differs from ρ_A(0) by less than 1.3·10⁻¹⁵ over t ∈ {1, 5, 10, 40}.
Multi-time correlations, responses to interventions and dynamics relative to a clock subsystem
are not frozen; these are exactly what linear-response physics in KMS states uses. One-time
local evolution needs a state rule that breaks equivariance.

### S10. What closure plus one hypothesis does force: equilibration on average

Theorem (Linden–Popescu–Short–Winter 2009): for a pure global state and an H with
non-degenerate energy gaps, ⟨D(ρ_A(t), ω_A)⟩_t ≤ ½√(d_A ⟨‖ρ_A − ω_A‖₂²⟩_t) ≤ ½√(d_A²/d_eff),
where ω is the time-averaged state and d_eff = 1/Tr ω². In this H the smallest level spacing is
5.9·10⁻⁴ and the closest pair of distinct Bohr gaps is separated by 3.2·10⁻⁸: the hypothesis
holds numerically, but the infinite-time average is reached only on times of order 10⁷, so the
sampled average over 3000 random times in [0, 10⁴] is an estimate.

| Pure preparation | d_eff | ½√(d_A²/d_eff) | ½√(d_A⟨‖·‖₂²⟩), exact average | sampled ⟨D⟩ | sampled max D |
|---|---|---|---|---|---|
| \|+⟩ ⊗ excited bath eigenstate (S8); d_eff < d_A² = 4, outside the non-trivial regime | 3.23 | 0.556 | 0.298 | 0.287 | 0.515 |
| \|+⟩ ⊗ bath Néel state \|0101010⟩ | 29.56 | 0.184 | 0.095 | 0.088 | 0.242 |

The bound holds, and its size is set by the supplied preparation through d_eff. Closure plus
generic gaps therefore forces equilibration **on average**, non-trivially when d_A² ≪ d_eff; it
does not force pointwise or monotone relaxation, a thermal ω_A, or a time scale.

### S11. Pointwise local inequalities that closure does force

**Athermality budget (closure alone, FORMAL_CORE THM 2.14).** With σ the Gibbs state of the
global H at β = 1 (stationary), D(ρ_Ω(t)‖σ) = 0.891092, constant to 4·10⁻¹⁵, and
D(ρ_A(t)‖Tr_Ā σ) ≤ 0.891092 at every one of the 41 unit-grid times (maximum 0.801675, at t = 0).

**Second-law identity (closure plus a Gibbs complement uncorrelated with A, COND 2.16).** The S1
preparation is ρ_A(0) ⊗ γ_Ā with γ_Ā the Gibbs state of H_Ā at β = 1, so
Σ(t) = ΔS_A + βΔE_Ā = I(A:Ā)(t) + D(ρ_Ā(t)‖γ_Ā) must hold exactly.

| Quantity | Value |
|---|---|
| max_t \|Σ − (I + D)\| over the unit grid | 1.1·10⁻¹⁴ |
| min_t Σ | −5·10⁻¹⁵ (zero at t = 0) |
| Σ at t = 1, 10, 20, 40 | 0.3457, 0.7391, 0.5513, 0.6648 |
| decreasing unit steps of Σ | 22 of 40 |

Σ is non-negative at every time and not monotone. The forced half of the arrow is this
pointwise inequality given a Gibbs reference preparation uncorrelated with A (general passivity is not sufficient); the direction is the
choice of that reference (S7).



## Calculation 2: nested and overlapping compatibility

| Check | Result |
|---|---|
| C1 nested restrictions of one global 4-qubit state: Tr_C(Tr_D ρ) and Tr_D(Tr_C ρ) vs Tr_{CD} ρ; ρ_B via AB vs via BC | identical in both trace orders (0; 3·10⁻¹⁷; 3·10⁻¹⁸) |
| C2 chain A–B–C, ρ_AB = ρ_BC = Φ⁺: overlap ρ_B = 1/2 from both sides | overlap agrees exactly; **no global state** (infeasible with CLARABEL and with SCS); smallest achievable Frobenius violation of the two marginals **0.2887** |
| C3 max F_AB + F_BC over all 3-qubit states (linear SDP) | **1.5000** (would be 2 if compatible); symmetric max 0.7500 each |
| C4 isotropic W(p) = pΦ⁺ + (1−p)1/4 on AB and BC: largest compatible p | **2/3** (0.666667); at p = 2/3 + 10⁻³ the minimum violation is 8.7·10⁻⁴ > 0 |
| C4 isotropic W(p) on AB, BC and CA (triangle): largest compatible p | **1/3** (0.333333), the separability threshold; violation 5.0·10⁻⁴ at 1/3 + 10⁻³ |
| C5 classical chain P(a,b), P(b,c): extension by Markov gluing | always exists (error 6·10⁻¹⁷) |
| C5 classical 3-cycle of perfect anticorrelations | no joint distribution (LP infeasible) |
| C5 classical 4-cycle of four perfect anticorrelations | joint exists (½ on 0101, ½ on 1010): a cycle permits failure but does not force it |
| C6 overlaps that disagree (ρ_B = diag(0.6,0.4) vs 1/2) | infeasible at once (both solvers); minimum violation 0.0500 |

Analytic proof for C2: a pure ρ_AB forces ρ_ABC = ρ_AB ⊗ ρ_C, hence ρ_BC = (1/2) ⊗ ρ_C, whose
fidelity with Φ⁺ is Tr(ρ_C)/4 = 1/4. The C3 optimum is the top eigenspace of
Φ⁺_AB ⊗ 1 + 1 ⊗ Φ⁺_BC (eigenvalue 3/2), e.g. (2|000⟩ + |110⟩ + |011⟩)/√6 with F_AB = F_BC = 3/4.
W(p) is the isotropic state; pairwise it is locally equivalent to the Werner state, and for the
chain globally (a Y on B), so 2/3 holds for both; the two triangle problems are distinct and
both give 1/3. Overlap agreement is therefore **necessary and not sufficient** for a common
global realization. Classically, overlap agreement suffices for every consistent family exactly
on acyclic access structures (Vorob'ev 1962); on a cycle a family fails exactly when it violates
an n-cycle inequality (Araújo et al. 2013), with an odd number of perfect anticorrelations as
the deterministic extreme;
quantumly it already fails on a chain (monogamy of entanglement). This is the quantum marginal
problem, standard since the 1960s; it is stated here because it is the exact form of the
working theory's requirement that every sector's effective environment be compatible with one
global model and preparation.

## Reading of the two calculations for the working theory

1. Closure alone forces the global theorems (unitarity, global conservation, entropy
   relations, quasi-recurrence), the restriction/compatibility rules, and some pointwise local
   inequalities: the entropy inequalities, monogamy and the athermality budget (S8, S11).
   Closure plus generic gaps forces equilibration on average (S10); closure plus a Gibbs,
   uncorrelated complement forces Σ = I + D ≥ 0 at every time (S11). All are standard.
2. Closure forces no local dynamical law: not monotone relaxation, not Markovianity, not a
   kernel, not a unique local law, not thermality of the complement, not the direction of the
   arrow, not record persistence. Each is exhibited by a legal closed model above. A unitarily
   equivariant state rule built from H alone is stationary, so its one-time marginals do not
   evolve (S9).
3. The reduced dynamics, kernel and noise are downstream of the supplied global state,
   Hamiltonian and access structure; identical marginals do not fix the future.
4. Records are subsystems; interventions are either externally timed (supplied controls, as
   in S6, which inject energy) or internal (a clock and controller inside Ω). Multi-time
   statistics are not functions of one-time maps, and a record the dynamics can write is not
   protected for every state (THM 2.15).
5. Compatibility of overlapping descriptions is a nontrivial constraint with a quantitative
   edge (3/2, 2/3, 1/3 above), and it is a standard constraint, not a GRUT discovery.

## Verification record

All records named below are in `results/verification/`.

- **Current version** (scripts and results as committed): `python3 code/verify.py --rerun`
  re-executes both calculations and requires identical leaf paths and agreement on every
  numeric leaf with the committed JSON (worst relative deviation 0 in this session). The
  separate agent's spec-only re-implementation was re-compared against the committed results
  after the final edits: 23 rows, all agree, script hashes matched
  (`independent_reimplementation_comparison_current.json`).
- **Independent re-implementation** (separate agent, written from the specification without
  reading the code; code `independent_reimplementation_code.py`, comparator
  `independent_reimplementation_compare.py`, series `independent_reimplementation_series.npz`):
  23 rows covering calculation-1 quantities (E_tot and its drift, the full 401-point E_tot(t)
  and S_A(t) series, S_A at t = 10 and 40 and its maximum, the backflow sum and D's minimum, the
  three record probabilities, the identical-marginal trace distance) and the calculation-2
  isotropic chain threshold. Twenty-one agree to 10⁻⁶ relative or better; the two numerical
  zeros (the E_tot drift and a zero trace distance) agree absolutely below 10⁻¹³; offsets are
  12-digit JSON rounding and solver tolerance. The comparison made after the first code review
  is kept as `independent_reimplementation_comparison_20oct_after_first_code_review.json`, and
  an earlier first attempt of the same agent as `independent_reimplementation_code_19oct_first_attempt.py`.
- **Independent re-execution** (separate agent, 19 Oct): verified the **pre-review** version of
  both scripts (hashes ba50f60e… and 1aad4b06…), 5030 and 24 numeric leaves bit-exact, two
  runs bitwise identical. That version is superseded; the record is kept as history
  (`independent_reexecution_19oct_prereview.md`).
- **Adversarial code review of calculation 1** (separate agent): "sound with minor issues", nine
  items, all addressed: entropy inequalities moved to the mixed trajectory; the S7 identity
  ρ_rev(t) = ρ(T−t)* recorded and the wording made non-monotone; the κX_R case relabelled a
  basis rotation and a bath-coupled register added; the composed one-time-map prediction
  recorded; the S6 identities labelled; the t = 0 qualifier on the inhomogeneity; the
  integrability claim withdrawn and the symmetries stated; the first local maximum quoted and
  the backflow sum labelled a grid quantity; dead code removed. The review confirmed the
  partial trace against brute force, the conditioning of the map inversion (at most 2.86 over
  the full grid), the step-robustness of the generator, the CP-divisibility failure at other τ
  and through negative canonical rates, the analytic inhomogeneity, the b′ period, the CNOT
  and the validity of the time reversal (`code_review_calculation1.md`).
- **Adversarial code review of calculation 2** (separate agent): "sound with minor issues",
  seven items (five prose defects, one statement confirmed accurate, one robustness caveat); the
  six actionable ones all addressed: the C2 fidelity
  corrected from 1/2 to 1/4; the C3 optimiser comment corrected; W(p) relabelled isotropic; the
  cycle note corrected; a second solver and minimum-violation certificates added; C1 checked in
  both trace orders. The review confirmed the cvxpy axis convention, the 3/2 and 3/4 optima by
  eigen-decomposition and 300 random restarts, the 2/3 chain threshold against the
  symmetric-extension criterion of Chen–Ji–Kribs–Lütkenhaus–Zeng (Phys. Rev. A 90, 032318
  (2014)), the 1/3 triangle threshold by a spin-decomposition construction, and the LP
  constructions (`code_review_calculation2.md`).
- **Attack on the packet**, round 1 (separate agents, equations / handoff compliance / number
  provenance): 2 blocking, 11 major and 22 minor findings, all addressed; they produced S4's
  fourth row, S9, S10 and the corrections above (`attack_round1/`). Round 2 (overclaim,
  steelman, equations, number provenance): no blocking, 30 major and 50 minor findings, all
  addressed; they produced S11, the S6 relabelling, the corrected S9 scope and the counts above
  (`attack_round2/`).
- Re-run of the M1 matched uniform-field control script (`uniform_field_control.py` in the M1
  review packet) reproduced `control_results.txt` line for line (8.6 s); see
  `results/m1_matched_control_rerun.txt`. No other old solver was reopened.

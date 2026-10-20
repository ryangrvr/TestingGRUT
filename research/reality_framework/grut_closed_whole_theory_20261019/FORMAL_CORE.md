# Formal core: the closed-whole framework

19 October 2026. Packet `grut_closed_whole_theory_20261019`. This file states the parent
question precisely, defines closure and the global-to-local rules, lists what those rules
force and what they leave free, and fixes the **theory datum** that every GRUT model instance
must supply. Every statement is labelled:

- **DEF** definition · **THM** theorem of standard quantum theory (reference given) ·
  **SUPPLIED** an input no rule here fixes · **NOT FORCED** with a countermodel ·
  **BY DEFINITION** settled by how DEF 1.1–1.8 are written, not by a countermodel ·
  **SCOPE** outside what the finite model can address · **COND** a theorem that holds only
  under stated extra hypotheses · **OPEN** not established.

Nothing in this file is claimed as new physics. The GRUT contribution is the organization:
every sector of the program is read as a campaign that asked the parent question in one
place, and the calculable framework below is what the campaigns' answers attach to.
Numbers quoted as "(calc 1, S n)" or "(calc 2, C n)" refer to [CALCULATIONS.md](CALCULATIONS.md).

## 0. The parent question

> What does a closed system, whose every subsystem's bath is the rest of the universe,
> force that known physics does not already contain?

An answer has the form: *in model class M, closure plus the stated supplied data implies X*,
together with either a derivation (X forced) or a countermodel (X not forced). "Nothing
additional was forced in class M, and this computation establishes it" is a complete answer.
"Closure" is defined in §1; "known physics" means quantum theory with a supplied dynamics,
state, access structure and measurement rule, plus whatever classical field theory a sector
imports.

## 1. Definitions

**DEF 1.1 (closed whole, finite case).** A triple (ℋ_Ω, H, ρ_Ω(0)) with ℋ_Ω finite
dimensional, H self-adjoint and time independent, ρ_Ω(0) ≥ 0, Tr ρ_Ω(0) = 1. Its law is
ρ_Ω(t) = U(t) ρ_Ω(0) U(t)†, U(t) = e^{−iHt/ħ}. There is no exterior: no reservoir, drive
or measurement device outside Ω. The time parameter t of U(t) is **SUPPLIED** and is read by
no subsystem of the model unless a clock subsystem is added (DEF 1.6).

**DEF 1.2 (access structure).** A specified family of subsystems. Finite case: a tensor
factorization ℋ_Ω = ℋ_A ⊗ ℋ_Ā for each subsystem A in the family, equivalently a family of
type-I factors 𝔐_A ⊂ B(ℋ_Ω). Algebraic case: a net O ↦ 𝔄(O) of von Neumann algebras on a
common Hilbert space with a global state ω. Which subsystems exist, and which operators an
observer can apply, is **SUPPLIED**. One relational route is known: a spectrum that admits a
k-local tensor-product description generically admits exactly one, up to local equivalences,
and most spectra admit none (Cotler–Penington–Ranard, Commun. Math. Phys. 368, 1267 (2019)).
If locality of H is supplied, it can therefore select the factorization; locality itself is
then the supplied input. No rule below derives the access structure from closure alone (the
program's EA0 campaign
asked whether it could be made endogenous; see the ledger).

**DEF 1.3 (bath = rest).** The complement of A is the commutant 𝔐_A′ (finite case: B(ℋ_Ā)).
This is a definition relative to the access structure. It carries no further property:
the complement need not be thermal, Markovian, infinite, weakly coupled, initially
uncorrelated with A, or resettable. Each of those is a **SUPPLIED** property of a model
instance when it holds (countermodels in §4).

**DEF 1.4 (restriction).** ρ_A = Tr_Ā ρ_Ω, the unique operator with Tr(ρ_A X) = Tr(ρ_Ω (X ⊗ 1))
for all X on ℋ_A. Algebraic case: ω_A = ω|_{𝔄(O)}. The restriction of a state exists for every
state and every subsystem; restricting does not select the state.

**DEF 1.5 (records and interventions).** A record is a subsystem R in the access structure
together with a unitary W on A ∪ R that correlates R with an observable of A. An intervention
at time t is an instrument {𝓔_x} on A (completely positive maps with Σ_x 𝓔_x trace
preserving). In a fully closed description an intervention is a unitary on A ∪ R with R an
ancilla in the access structure, followed by restriction. The complete multi-time object is
the process tensor 𝒯_{t_k…t_1} (Pollock et al., Phys. Rev. A 97, 012127 (2018); Phys. Rev.
Lett. 120, 040405 (2018)): a multilinear completely positive map from the interventions at
t_1 < … < t_k to the final state of A, determined by (H, ρ_Ω(0), access structure).

**DEF 1.6 (external time).** In the finite model t is a parameter of the one-parameter group
U(t); the clock that reads it is not part of Ω. Making the clock a subsystem (Page–Wootters,
Phys. Rev. D 27, 2885 (1983)) is an additional **SUPPLIED** construction with its own
assumptions; nothing here solves the problem of time. In general relativity there is no
preferred global t and, in general, no conserved total energy; the finite statements below
about energy conservation do not carry over without a stated asymptotic or gauge structure.

**DEF 1.7 (theory datum).** A model instance of the framework is the tuple

    Θ = ( 𝔄 or ℋ_Ω ,  H (or dynamics) ,  ρ_Ω(0) (preparation) ,  {𝔐_A} (access structure) ,
          {instruments, records} ,  clock / time parameter ).

A **prediction** is a Born probability Tr(ρ_Ω(t) E), E an effect in the access structure,
or a multi-time probability from the process tensor; every prediction is conditional on Θ.
Everything in Θ is supplied unless a rule below fixes it. (Calc 1 lists one complete Θ.)

**DEF 1.8 (measurement rule).** Outcome probabilities are given by the Born rule p(x) =
Tr(ρ_A E_x). This is **SUPPLIED** as a postulate (GRUT III P5), or traded for other
operational assumptions (Gleason 1957; Busch 2003; Masanes–Galley–Müller, Nat. Commun. 10,
1361 (2019), primary text inspected: their assumptions are the unitary-QM postulates for
states, transformations and composite systems, the "possibility of state estimation", and the
requirement that the set of outcome-probability functions be closed under mixtures, under
unitaries and under system composition (their Def. 9); these are operational assumptions about
measurements, not the dynamical closure of DEF 1.1). Dynamical closure does not derive it.

## 2. Exact global-to-local rules

**THM 2.1 (global invariants).** Under DEF 1.1 the spectrum of ρ_Ω(t), hence S(ρ_Ω) and
Tr ρ_Ω², and ⟨f(H)⟩ for every function f, are constant; any operator commuting with H is
conserved. (Unitarity. Calc 1, S1 and S8: drifts ≤ 10⁻¹⁴.)

**THM 2.2 (exact reduced equation).** For the decomposition H = H_A ⊗ 1 + 1 ⊗ H_Ā + H_{AĀ},

    dρ_A/dt = −(i/ħ)[H_A, ρ_A] − (i/ħ) Tr_Ā [H_{AĀ}, ρ_Ω(t)],

because Tr_Ā[1 ⊗ H_Ā, ρ_Ω] = 0 by cyclicity of the partial trace. The second term is a
functional of the **global** state, not of ρ_A. (Calc 1, S1: residual 1.9·10⁻¹⁰ by finite
differences; the H_Ā term 3·10⁻¹⁶.)

**THM 2.3 (Nakajima–Zwanzig and time-convolutionless forms).** With the projector
𝒫ρ = (Tr_Ā ρ) ⊗ ρ_ref, 𝒬 = 1 − 𝒫, and ℒ· = −(i/ħ)[H, ·],

    𝒫ρ̇(t) = 𝒫ℒ𝒫ρ(t) + ∫₀ᵗ ds 𝒦(t−s) 𝒫ρ(s) + 𝓘(t),
    𝒦(τ) = 𝒫ℒ e^{𝒬ℒτ} 𝒬ℒ𝒫,     𝓘(t) = 𝒫ℒ e^{𝒬ℒt} 𝒬ρ(0),

(Nakajima, Prog. Theor. Phys. 20, 948 (1958); Zwanzig, J. Chem. Phys. 33, 1338 (1960)), or
equivalently, wherever the relevant operator 1 − Σ(t) is invertible (for product
preparations, wherever Φ_t is invertible), ρ̇_A = 𝒦_TCL(t) ρ_A + 𝓘_TCL(t)
(Shibata–Takahashi–Hashitsume, J. Stat. Phys. 17, 171 (1977); Breuer–Petruccione, *The Theory
of Open Quantum Systems*, §9.2). Here ρ_ref is a density operator on ℋ_Ā, so 𝒫² = 𝒫.
Consequences: (i) the memory kernel depends on all of H, on ρ_ref and on the access structure;
(ii) the inhomogeneous term vanishes **if** 𝒬ρ(0) = 0, i.e. if ρ_Ω(0) = ρ_A(0) ⊗ ρ_ref; the
converse fails: a correlation χ with Tr_Ā χ = 0 and [H, χ] = 0 gives 𝒬ρ(0) ≠ 0 and 𝓘 ≡ 0,
while a product ρ_A ⊗ σ with σ ≠ ρ_ref has 𝒬ρ(0) ≠ 0 and generically 𝓘 ≠ 0; (iii) no
approximation (Born, Markov, secular, Gaussian) is implied. Kernel and noise are
**downstream** of the supplied datum. (Calc 1, S3: the exact TCL generator is time dependent
and differs between thermal, ground and infinite-temperature bath preparations; S4: the
A-factor of 𝓘(0) is ¼ Z₀ in closed form, and the conserved Z-parity correlation of the same
model gives identical futures, trace distance 4·10⁻¹⁶.)

**THM 2.4 (dynamical maps and their limits).** For product preparations with a fixed bath
state ρ_Ā, the map Φ_t: ρ_A(0) ↦ ρ_A(t) is linear and completely positive trace preserving
(Stinespring 1955; Kraus 1971). It need not be CP-divisible: there need not exist a
completely positive V_{t₂,t₁} with Φ_{t₂} = V_{t₂,t₁}Φ_{t₁} (Rivas–Huelga–Plenio, Phys. Rev. Lett. 105, 050403 (2010); Breuer–
Laine–Piilo, Phys. Rev. Lett. 103, 210401 (2009)). For correlated preparations there is in
general **no** map on ρ_A(0) that reproduces ρ_A(t) (Pechukas, Phys. Rev. Lett. 73, 1060
(1994); Štelmachovič–Bužek, Phys. Rev. A 64, 062106 (2001)): two global states with identical
ρ_A(0), and even identical ρ_Ā(0), can have different ρ_A(t). (Calc 1, S3: intermediate-map
Choi eigenvalue −0.159; S4: trace distance 0.1025 between the futures of two preparations
with identical marginals.)

**THM 2.5 (influence functional and the in-in generating functional are derived).** In the
path representation, Φ_t has the Feynman–Vernon form with influence functional
F[x, x′] = Tr_Ā[ U_Ā[x] ρ_Ā U_Ā[x′]† ] for a factorized preparation ρ_A ⊗ ρ_Ā (Feynman–Vernon,
Ann. Phys. 24, 118 (1963)); correlated preparations carry an additional preparation functional
and are not of this form (Hakim–Ambegaokar, Phys. Rev. A 32, 423 (1985); Grabert–Schramm–Ingold,
Phys. Rep. 168, 115 (1988)). **COND:** if the complement is linear (harmonic), coupled
linearly, and prepared in a factorized Gaussian state, F is Gaussian and is fixed by two kernels built from a
spectral density J(ω): the dissipation kernel η(τ) ∝ ∫ dω J(ω) sin ωτ and the noise kernel
ν(τ) ∝ ∫ dω J(ω) coth(βħω/2) cos ωτ for a thermal preparation (Caldeira–Leggett, Physica A
121, 587 (1983)). The closed-time-path / in-in generating functional with doubled fields
(Schwinger 1961; Keldysh 1964) is the same object in field-theoretic notation. Therefore a
quadratic in-in action S_CTP is the Gaussian approximation to a derived functional; it
presupposes J(ω), β, linearity and the preparation. It is not a parent principle. This places
the parent object that the handoff attributes to GRUT V7 (the V7 file itself is not in the
repository): downstream, conditional, and standard.

**THM 2.6 (nested compatibility).** If A ⊂ B then Tr_{B∖A} ρ_B = ρ_A for the restrictions of
one global state; restriction maps compose. Automatic. (Calc 2, C1.)

**THM 2.7 (overlap compatibility is necessary and not sufficient).** For subsystems A, B with
A ∩ B ≠ ∅, restrictions of one global state agree on A ∩ B (necessary). Conversely, a family
{ρ_{A_i}} that agrees on all overlaps need not be the restriction of any positive global
state: with ρ_AB = ρ_BC = |Φ⁺⟩⟨Φ⁺|, the overlaps agree (ρ_B = 1/2 from both) and no ρ_ABC
exists, because a pure ρ_AB forces ρ_ABC = ρ_AB ⊗ ρ_C and then ρ_BC = (1/2) ⊗ ρ_C, whose
fidelity with Φ⁺ is 1/4 (monogamy of entanglement; the quantitative tangle inequality is
Coffman–Kundu–Wootters, Phys. Rev. A 61, 052306 (2000), not needed here). Existence of a
common extension is the quantum marginal problem (a semidefinite feasibility problem in
finite dimension; QMA-complete in general under Turing reductions, Liu 2006, and under Karp
reductions, Broadbent–Grilo, SIAM J. Comput. 51, 1400 (2022); Klyachko 2006 for the one-body
pure case). Classically, overlap agreement suffices for every consistent family iff the
hypergraph of subsystems is acyclic (Vorob'ev, Theory Probab. Appl. 7, 147 (1962)); on a
cycle a consistent family fails exactly when it violates an n-cycle inequality (Araújo et
al., Phys. Rev. A 88, 022118 (2013)), an odd number of perfect anticorrelations being the
deterministic extreme, and need not fail otherwise; quantumly failure occurs already on a
chain.
(Calc 2: Bell chain infeasible, minimum marginal violation 0.289; max F_AB + F_BC = 3/2;
isotropic-state chain threshold 2/3, triangle 1/3; classical chain always extends, classical
odd anticorrelated 3-cycle fails, 4-cycle extends.)

**THM 2.8 (records and multi-time statistics).** A record written unitarily into an idle
register R and then traced out yields, for all later statistics of A, the same result as the
corresponding instrument applied to A; and the statistics after an intervention are not
functions of the one-time map Φ_t (back-action). (Calc 1, S6: Born probability 0.354910 =
register population; p(Z=+1, t₂) 0.3385 with the record vs 0.4491 without.) The multi-time
object is the process tensor of DEF 1.5; an operational Markov condition is a property of it,
not of Φ_t (Pollock et al. 2018).

**THM 2.9 (entropy relations and recurrence).** For a pure ρ_Ω: S(ρ_A) = S(ρ_Ā).
In general: subadditivity, strong subadditivity (Lieb–Ruskai 1973), Araki–Lieb
|S_A − S_Ā| ≤ S_Ω ≤ S_A + S_Ā. For finite dimension, t ↦ ρ_Ω(t) is almost periodic, so the
state returns arbitrarily close to ρ_Ω(0) (quasi-recurrence; Bocchieri–Loinger, Phys. Rev.
107, 337 (1957)); no bound on the recurrence time is implied. (Calc 1, S8: identities on a
pure global state; subadditivity and Araki–Lieb checked non-trivially on the mixed S1
trajectory, where they are strict after t = 0.)

The next four statements are **conditional**: closure plus one named, checkable hypothesis
forces a local consequence. They are standard theorems; the hypothesis in each is a property
of the supplied datum, not of closure.

**COND 2.10 (equilibration on average).** If H has non-degenerate energy gaps (E_k − E_l =
E_m − E_n only trivially) and ρ_Ω is pure, the infinite-time average of the distance of ρ_A(t)
from its time average ω_A obeys ⟨D(ρ_A(t), ω_A)⟩_t ≤ ½ √(d_A² / d_eff), with d_eff = 1/Tr ω²
the effective number of energy eigenstates in the preparation (Linden–Popescu–Short–Winter,
Phys. Rev. E 79, 061103 (2009), Theorem 1; Reimann, Phys. Rev. Lett. 101, 190403 (2008);
mixed states and degenerate gaps with modified constants: Short, New J. Phys. 13, 053009
(2011); Short–Farrelly, New J. Phys. 14, 013063 (2012)). This forces equilibration **on
average**, for every subsystem small compared with d_eff. It does not force pointwise or
monotone relaxation, a thermal ω_A, or any time scale, and its size is set by the supplied
preparation through d_eff. (Calc 1, S10: the closest pair of distinct Bohr gaps is separated
by 3·10⁻⁸, so the hypothesis holds numerically but the average converges only on times of
order 10⁷; for the S8 preparation d_eff = 3.2, bound 0.556, sampled ⟨D⟩ = 0.287; for a Néel
bath preparation d_eff = 29.6, bound 0.184, sampled ⟨D⟩ = 0.088.)

**COND 2.11 (canonical typicality).** If the global state is drawn from the unitarily
invariant measure on a subspace ℋ_R of large dimension, then with overwhelming probability
ρ_A is close to Tr_Ā(1_R/d_R) (Popescu–Short–Winter, Nat. Phys. 2, 754 (2006);
Goldstein–Lebowitz–Tumulka–Zanghì, Phys. Rev. Lett. 96, 050403 (2006)). The measure over
preparations is the supplied hypothesis; closure does not choose it.

**COND 2.12 (effective light cones).** If H is a sum of bounded local terms on a lattice, local
operators spread with a finite Lieb–Robinson velocity up to exponentially small tails
(Lieb–Robinson, Commun. Math. Phys. 28, 251 (1972)). Locality of H is the supplied hypothesis.

**COND 2.13 (conservation limits on records).** If H conserves an additive quantity L, an
observable that does not commute with L cannot be measured exactly and repeatably by a unitary
coupling to an apparatus inside the whole; the error falls as the apparatus's spread in L
grows (Wigner 1952; Araki–Yanase, Phys. Rev. 120, 622 (1960); Ozawa, Phys. Rev. Lett. 88,
050402 (2002)). Here closure matters directly, because every apparatus is inside Ω; the
conserved quantity is a property of the supplied H.

## 3. What closure forces

From closure alone: unitarity and its invariants (2.1), the exact reduced equation (2.2), the
NZ/TCL structure with its inhomogeneity (2.3), CPTP maps for product preparations and their
possible non-divisibility (2.4), the derived status of influence functionals and of S_CTP
(2.5), nested compatibility (2.6), necessary overlap compatibility (2.7), the record/
instrument equivalence and the multi-time structure (2.8), the entropy relations and
quasi-recurrence (2.9), no-signalling between commuting subsystems, and monogamy constraints.

From closure plus one named hypothesis on the supplied datum: equilibration on average for
non-degenerate gaps (2.10), canonical typicality for a stated measure over preparations
(2.11), effective light cones for local H (2.12), and conservation-law limits on records made
inside the whole (2.13).

All are standard. None is a GRUT contribution; the GRUT content is that each sector is held
to exactly these rules and nothing more unless it supplies more.

## 4. What closure does not force (each with a countermodel or a stated reason)

| Claim "closure forces …" | Status | Countermodel or reason |
|---|---|---|
| a pure global state | BY DEFINITION | DEF 1.1 admits mixed ρ_Ω, and calc 1 S1 uses one with all rules holding. An argument that a whole with no exterior must be pure is excluded by this definition rather than refuted by a countermodel |
| a preferred initial state | NOT FORCED | any ρ_Ω(0) is admissible; calc 1 S7: the time-reversed final state is a legal preparation; calc 1 S9: rules F(H) that depend on H alone freeze all local dynamics |
| a thermal complement | NOT FORCED | calc 1 S3: ground, thermal and infinite-temperature complements all legal, with different kernels |
| pointwise or monotone relaxation | NOT FORCED | calc 1 S2 (26 sign changes of dS_A/dt, backflow 1.20), S5 (exact return at T = 4π). Equilibration **on average** is forced under non-degenerate gaps (COND 2.10), with size set by d_eff |
| a Markovian or exponential memory kernel | NOT FORCED | calc 1 S3: exact TCL generator time dependent, map not CP-divisible |
| a unique local law for A | NOT FORCED | calc 1 S5: same H_A and ρ_A(0), three complements, trace distances up to 0.79 |
| that ρ_A(0) determines ρ_A(t) | NOT FORCED | THM 2.4; calc 1 S4 |
| persistence of records | NOT FORCED | records are stable for all states iff [H, Z_R] = 0; calc 1 S6: otherwise retrievability depends on the supplied coupling (lossless for a register rotation, reduced for a bath-coupled register) |
| an arrow of time | NOT FORCED | calc 1 S7: S_A returns from 0.42 to 0 along a legal preparation; the arrow is the low-entropy preparation (GRUT III P4 locates it as a postulate) |
| a small gravitational Lyapunov exponent or any chaos exponent | SCOPE | in finite dimension t ↦ ρ_Ω(t) is almost periodic (THM 2.9), so no exponential growth persists for all times; out-of-time-order growth at intermediate times is model specific; no gravitational sector is in any model here |
| compact spatial topology, finite four-volume, positive spatial curvature | SCOPE | "closed system" (no exterior) is not "closed universe" (k = +1); spatially open FLRW models have no exterior either; the finite model has no geometry at all |
| the observed constants or a specified matter/gravity action | BY DEFINITION | they are elements of Θ (GRUT III data and postulates P0, P2, P3) |
| global energy conservation in cosmology | SCOPE | THM 2.1 needs a time-independent H; an expanding spacetime has none in general without asymptotic structure |
| the Born rule | BY DEFINITION | DEF 1.8 supplies it; dynamical closure is not one of the assumptions that the reconstruction theorems replace it with |

## 5. Where the finite picture stops

**COND 5.1 (algebraic setting).** In quantum field theory the local algebras 𝔄(O) are of type
III, and of type III₁ under a scaling-limit hypothesis (Fredenhagen, Commun. Math. Phys. 97,
79 (1985); Buchholz–D'Antoni–Fredenhagen, Commun. Math. Phys. 111, 123 (1987)); there is no
tensor factorization ℋ = ℋ_O ⊗ ℋ_{O′}, no partial trace and no reduced density matrix (review:
Witten, Rev. Mod. Phys. 90, 045003 (2018), §§6–7). DEF 1.3–1.4 must be read
algebraically: bath = commutant 𝔄(O)′, equal to the causal-complement algebra 𝔄(O′) only
where Haag duality holds; restriction = ω|_{𝔄(O)}. GRUT III's P1 is exactly this
identification and supplies no dynamics and no state (GRUT III, 01_PRINCIPLES; closure-fork
response). A regulated subregion density matrix requires an explicit split or cutoff.

**COND 5.2 (modular flow).** For a faithful normal state the restriction is KMS with respect
to its own modular flow (Tomita–Takesaki): at β = −1 for σ_t(A) = Δ^{it}AΔ^{−it}, equivalently
at β = 1 for σ_{−t} (Bratteli–Robinson II, Thm 5.3.10; sign conventions differ between
authors). This is automatic and selects nothing;
it is a temperature in a physical clock only where the modular flow is geometric
(Bisognano–Wichmann; Kay–Wald; Gibbons–Hawking), which are imports. Equating modular time
with physical time in general is an additional hypothesis (Connes–Rovelli, Class. Quantum
Grav. 11, 2899 (1994)), not a consequence of closure.

**COND 5.3 (natural states).** Primary text inspected: Fewster–Verch, Ann. Henri Poincaré
13, 1613 (2012), Definition 6.12 and Theorem 6.13. A natural state is an assignment M ↦ ω_M
with A(ψ)*ω_N = ω_M for **every** morphism ψ: M → N of the spacetime category Loc. Theorem
6.13: if the theory is dynamically local and has a natural state, and if in some spacetime M
with **noncompact** Cauchy surfaces ω_M induces a **faithful** GNS representation with the
**Reeh–Schlieder** property, then the relative Cauchy evolution in M is trivial; with extended
locality in addition, the theory is equivalent to the trivial theory. The noncompactness is
used in the proof. The paper describes the result as bringing to sharper form an argument
sketched by Brunetti–Fredenhagen–Verch (2003) for the free scalar field. Scope: the theorem
says nothing about rules covariant only under isometries, about finite models, or about a
category of spatially compact spacetimes; it excludes one specific kind of covariant
state-selection rule under these hypotheses. In the finite class the corresponding
obstruction is elementary: a rule F(H) equivariant under unitary conjugation commutes with H,
so the selected state is stationary and every subsystem is frozen (calc 1 S9). A
state-selection rule with local dynamical content must therefore depend on more than H, and
the program's C1/C2 and CP-12 campaigns are the record of attempts.

**COND 5.4 (gravity).** With diffeomorphism constraints there are no local gauge-invariant
subsystems without dressing, no preferred t, and in general no conserved total energy.
Nothing in this packet addresses these; the finite statements are the regulated model's.

## 6. The theory datum as a checklist

Before any prediction is stated, a GRUT model instance must list, in this order:

1. the algebra or Hilbert space, with its regulator if continuum;
2. the dynamics: H, or an action with its field content and couplings;
3. the preparation ρ_Ω(0) (or ω), and whether it is pure, product, Gaussian, thermal;
4. the access structure: subsystems, their factorization or algebras, the observables
   available to each observer;
5. the instruments and records, including registers that store outcomes;
6. the clock or time parameter and its relation to the dynamics.

A sector that states a "response", "kernel", "temperature", "relaxation rate" or
"prediction" without items 2–4 has not yet stated a prediction. A sector that supplies all
six and computes a Born probability has stated a prediction **conditional on Θ**.

## 7. What an additional law would have to be

A GRUT postulate with content beyond §2 must constrain Θ. The slots of Θ (DEF 1.7) are
exhaustive by construction; a law may constrain one slot, relate two slots, impose a measure
over Θ, or modify the rules of DEF 1.1 and DEF 1.8 themselves. The table lists each kind.

| Type | What it would fix | Cost | Program campaigns that tested it |
|---|---|---|---|
| algebra slot | the Hilbert space, field content, gauge algebra or regulator | must recover the observed spectrum; GRUT III treats it as data | GRUT III P0 and data; CP-12 graph Hilbert space; Candidate 1 field content |
| state selection | ρ_Ω(0) = F(H, access, …) | in the finite class a rule equivariant in H alone is stationary (calc 1 S9), so a rule with content must use the access structure, a clock or a non-equivariant datum; in locally covariant QFT an embedding-covariant rule is excluded under the hypotheses of COND 5.3; must give a positive normalizable state | closure / no-boundary (C1, C2, B control), CP-12 preparation, GRUT I–II projector choice |
| dynamical constraint | restrictions on H or the action (e.g. universal coupling) | must recover known physics; introduces its own constants | GRUT III P2–P3, Candidate 1, influence hierarchy |
| access constraint | which subsystems and observables exist, endogenously | must not presuppose the factorization it derives | EA0 endogenous access, Level-0 floor, CP-12 relational record gates |
| relational constraint between slots | e.g. locality of H selecting the factorization (Cotler–Penington–Ranard 2019) | the relation itself is supplied (locality); generic H admit no local factorization | EA0 and the Level-0 access bridge (NONUNIQUE-LIFT) |
| instrument / record slot | which interventions and records exist inside the whole | conservation limits (COND 2.13); records are subsystems whose persistence is set by H (calc 1 S6) | M1 writer, B2 receiver, CP-12 record gates |
| clock slot | modular time = physical time, or a clock subsystem | COND 5.2; must reproduce geometric cases | GRUT III P1, two-temperatures study |
| measure over Θ | typicality or a past-hypothesis measure over preparations | COND 2.11 shows such a measure has consequences; the measure is the postulate | CP-12 min-entropy obligation; GRUT I–II arrow (direction imported) |
| modified global law | non-unitary or stochastic ρ̇_Ω | changes THM 2.1; empirically bounded | not a GRUT campaign (collapse models) |
| measurement rule | derive DEF 1.8 | must replace it by other assumptions | GRUT III P5 (withdrawn HJW derivation) |

A bottom-up proposal that assigns local states sector by sector must also pass the marginal
problem (THM 2.7); a top-down global state passes it automatically.

Each row's result is reported in the campaign ledger. The baseline against which any
addition is scored is this file: predictions conditional on Θ, with Θ fully supplied.

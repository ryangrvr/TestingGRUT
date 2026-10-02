# SCOUT-2 ZOOM-OUT 04 (after S2-Σ + G1)

**Central question (owner):** is Σ genuinely selected, or have locality, classicality and objectivity merely become a
multi-objective preference problem?

## Answer

**A multi-objective preference problem, plus two exact degeneracies.**

- **No TPS dominates without a scoring convention.** This held in 58/58 runs over 4 models, 7 states, and a grid of
  τ, δ, fragment definitions and ε.
- The only scale-free object is the Pareto front, and it holds 3 – 148 inequivalent TPSs spanning 5 – 10 local
  dimension types.
- Every scalar winner depends on weights, normalization, transform, scale and state.

**Selection appears only by supplying Σ's own ingredients:** local dimension, number of factors, and a locality
framework. Even then it is up to the commutant of H. The state removes the commutant only by fixing an epoch.

**Verdict: Σ SURVIVES AS A SUPPLIED OBJECT** (a residual specification layer). Its irreducible content in these toy
models is:

> **Σ = (local dimension d, number of factors, objective family / preference order, scale (τ, δ, fragments), epoch)**
> **modulo the gauge of H-symmetries and the commutant of H.**

## The ten questions

1. **Derived or renamed?** Nothing derived. The literature conflict was **resolved by scope**:
   - CPR: unique given (k, d), up to the commutant;
   - Carroll–Singh: needs ψ, τ and a criterion;
   - Stoica: H-only selection is non-unique, shown exactly via the commutant.
2. **Where did the measure enter?** The Dirichlet weight sampling, the GUE perturbation, and the finite candidate
   family. All are named; none is earned.
3. **Where did access enter?** Fragment definition (R), coarse threshold δ, and the autonomy notion P. These are
   A_res-type items now located **inside** Σ's scale component.
4. **Is the subsystem decomposition unique?** No. In CPR scope it is unique only up to the commutant.
5. **Is composition unique?** Unchanged.
6. **Did a basin sneak in?** Yes: ψ and the epoch.
   - The state picks the winner.
   - In eigenstates and Gibbs states the state criteria go blind (Q ≡ 0), and selection falls back on H-only
     criteria, which are commutant-blind.
7. **Is an observer assumed?** Implicitly, through the preference order (which criterion matters) and the scale.
8. **Did information decrease?** No. Σ gained explicit components (epoch, preference order) and lost none.
9. **Has the C5 boundary moved?** Σ has been **characterized** rather than moved. It is now a specific tuple with a
   specific gauge group.
10. **Hostile that would overturn this:** a criterion family that is (i) commutant-sensitive without a state, or
    (ii) state-based but epoch-covariant in a way that selects a *time-independent* TPS, and (iii) dominant over all
    other admissible criteria without weights. Candidate for a future probe: trajectory-integrated criteria (averaged
    over all epochs), which remove the epoch dependence but reintroduce a time measure (H / measure-priced).

## Decomposition status

    C5 → D_dyn ⊕ H ⊕ Σ ⊕ A_res       (T2-2, classification; unchanged in form)

    Σ = (d, #factors, preference order, scale, epoch) / (Sym(H) × commutant(H))

    H splitting: H_measure → D (S2-3b, rigid compact D); H_preparation → A (S2-3b);
                 H_state now also carries the epoch / ψ dependence of Σ.

**Candidate T2-3** (classification; **CANONICAL-CANDIDATE: NO**): *H-only subsystem selectors are blind to the
commutant of H. Any selector that breaks this blindness uses the state at a chosen epoch.* This is exact in the tested
models (Σ-7b, numerical difference ≤ 4·10⁻¹⁶). It is a general identity for frames commuting with H, not merely a
toy-model observation. The general part is KNOWN-RESULT-adjacent (Stoica-type). The SCOUT-2 content is its
quantitative role in the Pareto analysis.

## Next

| Probe | Why |
|---|---|
| **S2-H2** (now unlocked) | attack H independently: a unique global attractor from every admissible initial condition with no supplied measure; contractive / primitive-Markov / gradient vs unitary; firewall: invertible dynamics cannot erase fine-grained information |
| S2-Σb (optional) | epoch-integrated criteria: does averaging over all epochs remove the epoch at the price of a time measure? |
| S2-G2 / G3 / G4 | graph-spectral, spacetime and capacity dimension, kept separate |

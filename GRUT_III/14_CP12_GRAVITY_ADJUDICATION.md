# 14 — Adjudication: CP-12 gravity response, the pair-tensor audit, and the proposed relational tensor charter

**Status (2026-10-13).** Adjudicates three documents: the CP-12 gravity-response report (`codex/grut-cp12-gravity-response`,
commit `11fbaec`), the same-day pair-tensor symmetry audit (`codex/grut-cp12-pair-tensor-symmetry-audit`, `a6d3823`), and
ChatGPT's review proposing a "CP-12 Relational Tensor Propagation and Universality Audit". What was checked:
- the report's `verify.py` was run on a full extract of the commit: **ten groups pass**. It needs scipy, which the manifest
  does not list, and it fails under `python3 -I` because it imports sibling modules;
- `exact_short_time.py` was run: **exact rational PASS** on path, cycle, star, ring-6 and ring-7;
- an independent re-implementation of the cubic response (different fermion ordering, own state builder), an orbit
  enumeration, a Γ-scan of the pair shell, and probe-overlap diagnostics (two verification agents);
- the literature added below, at full text where it matters (citation corrections are listed at the end).
No repository law or result was changed; this file is additive.

## 1. Verdicts on the report and the review

| Claim | Verdict |
|---|---|
| Absent-edge response R_ei(T) = Γ²t²√d_i T³/M² + O(T⁵), independent of distance | **Confirmed independently.** The T⁰, T¹, T², T⁴ terms vanish; occupied edges not incident on the source start at T⁵; the sum rule over all edges holds |
| Incident occupied edge: −Γ²t²(M−Q)T³/(M²√d_i) | **Confirmed** |
| Coefficient falls as N⁻⁴ for bounded source degree | **Confirmed** (M² ≈ N⁴/4) |
| Prepared geometry and ground state live in different harmonic sectors | Consistent with the record; not re-derived here |
| Pair-shell multiplicity N(N−3)/2 is an event-permutation property, not helicity | **Confirmed** |
| ChatGPT's theorem R_ij(T) = R_ik(T) for j, k ≠ i in a permutation-invariant state | **Correct.** It is the same result as the codex audit's three-orbit theorem, derived independently the same day |
| Orbit counts 1, 2, 4 with 0, 1, 2 events fixed | **Confirmed** for N = 6, 7, 8; in general C(f,2) + f + 1 orbits (needs N − f ≥ 2) |
| Pair gaps rise slightly over N = 5–7 and are not a massless mode | **Confirmed** |
| Favourable results come from different Q sectors | **Confirmed in detail**: the Q = 2N ground-state theorem is negative (not concentrated on 4-regular graphs); the pair shells are Q = N − 1 ground states with no locality result; the operational window holds only for prepared states with supplied geometry. No two favourable properties are shown in the same sector and the same state |

## 2. Two things the report and the review both understate

**The connection sector cannot create geometry locally.** The rewiring Laplacian is the only term that changes which
connections exist, and its off-diagonal amplitude, −Γ/M, is the same for every removable–absent pair whatever their graph
distance. Creation of a connection is therefore blind to matter. Matter acts on geometry only by *gating removal* (a
charged or flux-carrying connection cannot be removed: the T²-versus-T⁴ locking) and through the rewiring diagonal, which
counts removable edges. That is why the absent-edge response is distance-independent at every order where it is computed:
it is a global redistribution of a conserved resource, not a field sourced at a point. A Newtonian or gravitational
response that falls off with distance cannot come from this channel. If it exists in CP-12, it must come from matter
propagating on a prepared graph and then locking connections, so its locality is inherited from the prepared graph.

**The low-lying "tensor" shell is a slow connection mode, not a tensor response.** The two objects the report places side
by side are different:
- the *electric-pair tensor responses* are matter–electric excitations of a fixed, supplied graph, with gaps of about 2.5–3.1;
- the *pair shell* (gap ≈ 0.07 at N = 5–7) is the first Johnson harmonic of the connection sector, dressed by matter. A Γ
  scan at N = 5, Q = 4 gives gaps 0.0118, 0.0294, 0.0685, 0.149 at Γ = 0.025, 0.05, 0.1, 0.2: the gap scales with the rewiring
  rate and vanishes with it. An edge-occupation probe in the pair representation puts 99.5% of its variance in the shell.
  The electric-pair probe of the report's §7 couples to it with weight 3.9 × 10⁻⁹, and the T2 probe with 0.54%.
So the low shell is a nearly flat, slow rearrangement mode of the all-pairs connection pool, with no wave vector. Its
smallness is the smallness of Γ, not the approach to a massless pole. The report's phrase "their lowest common fivefold
shell" should be tightened so a reader does not take the shell for the tensor channel.

## 3. The charter: right structure, four constraints it must carry

ChatGPT's three questions (relational distance, tensor propagation, universal coupling) are the right ones, and its
same-regime requirement is essential. Its target is
χᴿ_ij,kl(ω, k) ∼ Z Π^TT_ij,kl(k) / (c_g²|k|² − (ω + i0)²). That target walks into four established results, which the charter
should state up front so a success cannot be claimed past them.

1. **Weinberg–Witten.** In a Poincaré-invariant theory with a conserved, Lorentz-covariant stress tensor that generates
   translations, no massless particle of helicity |h| > 1 can exist, elementary or composite (Phys. Lett. B 96 (1980) 59).
   CP-12 itself is not Lorentz invariant (external time, no boosts), so the theorem does not bite microscopically. It bites
   exactly on the regime the charter is aiming for: a Lorentzian infrared with a massless TT pole that couples universally to
   a conserved T. The known escape routes (Jenkins, IJMPD 18 (2009) 2249) are: (a) no stress tensor at all; (b) a fundamental
   gauge invariance under which T is not invariant, as in GR; (c) an emergent spacetime distinct from the background, as in
   AdS/CFT; (d) gravitons that stay non-relativistic. Jenkins argues that an emergent spin-2 gauge invariance in the
   background spacetime alone is not a consistent route. CP-12 is in (d) by default. **The charter must name its escape route
   before computing the pole.**
2. **Marolf (PRL 114, 031104, 2015).** Nonlinear gravity with a Gauss law for energy cannot emerge from a kinematically local
   theory, on a lattice or in the continuum, with or without Lorentz or translation invariance; strictly linear spin-2 modes
   are not constrained. The report's §8 already shows CP-12's U(1) Gauss law is not a Gauss law for energy. So question 3's
   "gravitational constraint emerging" needs a stated source of kinematic nonlocality for energy. Until then the realistic
   ceiling is a linear spin-2 mode.
3. **GW170817.** The speeds of gravity and light agree to −3 × 10⁻¹⁵ ≤ (v_GW − v_EM)/v_EM ≤ +7 × 10⁻¹⁶ (Abbott et al., ApJL 848
   (2017) L13; with stated emission-time assumptions). The charter's c_g must equal the matter light-cone speed to that
   precision. In CP-12 the connection sector does not propagate at all (§2), and matter moves with t_hop. The equality has to
   come out of the mechanism, not be tuned.
4. **Lattice gravitons already exist, and stop exactly where the charter's target begins.** Controlled lattice and qubit
   models give gapless helicity-±2 modes with ω ∝ k² (Xu, PRB 74, 224433 (2006): two tensor modes plus a scalar; Xu–Hořava,
   PRD 81, 104033 (2010), z = 2) or ω ∝ k³ (Gu–Wen L-type, Nucl. Phys. B 863 (2012) 90; Xu–Hořava z = 3, with only the two TT
   modes). Linear dispersion appears only in an uncontrolled approximation (Gu–Wen N-type, which its authors later called
   unreliable) or with a term added by hand (Xu). Xu–Hořava show the lattice symmetries of their model forbid z = 1, the
   Einstein–Hilbert term and the nonlinear self-interaction. None of these models has a nonlinear completion or universal
   matter coupling. **A linear-dispersion TT pole is the step no lattice model has reached in a controlled way**, so the
   charter should measure the dynamical exponent z and pre-register z = 1 as success; z ≥ 2 places CP-12 in the existing
   lattice-graviton class, which is not general relativity.

On question 1: the permutation theorem says a permutation-invariant state cannot encode distance in labelled two-point
responses. That is the discrete counterpart of a standard lesson: coordinate-labelled local quantities are not gauge
invariant, and only relational or label-averaged quantities are observable (Giddings–Marolf–Hartle, "Observables in
effective gravity", PRD 74 (2006) 064018; in causal-set theory, the label-invariant "stem" observables of Brightwell, Dowker,
García, Henson & Sorkin 2003). It is not a discrete analogue of Torre's theorem (no local-functional observables in closed
vacuum GR), and causal-set label invariance uses natural labellings only, a weaker symmetry than full permutations. The
practical consequence is the one both reviews draw: distance must be read from correlations conditioned on physical
records, or from a phase in which permutation symmetry is broken.

## 4. Recommended amendments (paste-ready)

```
CP-12 Relational Tensor Propagation and Universality Audit — amendments (Claude)
0. One sector, one state: declare Q/N and the state (prepared / ground / mixture)
   once; all three questions must be answered there. (Today no two favourable
   CP-12 results share a sector and a state.)
1. Relational distance: in an S_N-invariant state labelled two-point responses
   are uniform (C(f,2)+f+1 orbits with f fixed events). Use record-conditioned
   shells, V(r), spectral dimension; or exhibit S_N breaking. Note: connection
   creation is matter-blind (uniform -Gamma/M); matter only gates removal, so any
   distance dependence is inherited from a prepared graph.
2. Tensor propagation: separate the slow connection harmonic (gap ~ Gamma, 99.5%
   occupation-probe weight, 4e-9 electric-probe weight) from electric tensor
   responses (gaps ~2.5-3). Measure the dynamical exponent z; success is z = 1
   with exactly two TT channels. z >= 2 = existing lattice-graviton class
   (Xu 2006; Xu-Horava 2010; Gu-Wen 2012), not GR.
3. Universal coupling: name the Weinberg-Witten escape route (no T / fundamental
   gauge invariance / emergent spacetime / non-relativistic graviton) and the
   source of kinematic nonlocality for energy (Marolf 2015) BEFORE computing the
   pole. Compare c_g with the matter light-cone speed against GW170817
   (-3e-15 <= dv/v <= 7e-16); no tuning.
4. Keep geometry selection running in parallel; a tensor pole on an imposed
   graph does not explain why spacetime exists.
Reproducibility: verify.py passes 10/10 with scipy installed (not in MANIFEST)
and without python -I (sibling imports).
```

## 5. Citation corrections for the next writer

- Loebbert's Weinberg–Witten essay is Annalen Phys. 17 (2008) 803 and has **no arXiv number**; arXiv:0804.2193 is an
  unrelated paper. Open reviews: Jenkins, hep-th/0607239 and arXiv:0904.0453.
- Xu, PRB 74, 224433 (2006) is arXiv:cond-mat/0609595; cond-mat/0602443 is a different, shorter paper.
- Torre, PRD 48 (1993) R2373 is titled "Gravitational observables and local symmetries"; Giddings–Marolf–Hartle, PRD 74 (2006)
  064018 is titled "Observables in effective gravity".
- Gu & Wen, gr-qc/0606100 was never published in a journal; its linear-dispersion model is the "N-type" model their 2012
  paper calls unreliable.

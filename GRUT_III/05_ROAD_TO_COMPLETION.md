# 05 — The road to completion (revision 2): what is open, posed sharply, with the next computation

GRUT III has five framework postulates, a list of data, and — after the reviews — the admission that in its own
(KMS) regime it is standard physics. A theory of reality would derive the postulates from fewer and the data from
none, and would contain at least one statement conventional theory leaves free. Here is the distance, item by item.

## O1 — The area-law coefficient: where G comes from (central; made precise in `v9`, D12)

**Status.** P2 inserts η = 1/4G. In the heat-kernel scheme the species-by-species matching δ(1/4G) = (S_ent +
S_contact)/A is a theorem for free matter loops (minimal scalar 1/48πε², Dirac 1/24πε², conformal scalar 0 after the
Wald term, Maxwell −1/12πε² after Kabat's term). The coefficients are scheme-dependent; the scheme-free statement is
S_ent = ¼δ(M_P²)A (Cooperman–Luty). **With Standard-Model content the one-loop induced 1/G is +1/6 (ξ = 0) or
−1/2 (ξ = 1/6) in units of Λ²/π: tiny or negative.** Sakharov's route does not give G from this field content; the
finite G needs either a bare term, beyond-Standard-Model content (positive-k₁ species), or a UV completion in which
S = A/4G is exact (as Ryu–Takayanagi in AdS).
**Next computation.** None at one loop; the route is closed at this order. What remains is a demand on the UV
completion: exact horizon thermality, an area law with finite coefficient, P0 below the cutoff.

## O2 — The value of Λ

Integration constant (D3). Dead end recorded: the static-patch identity F = E − TS = 0 holds for every H and selects
nothing. Position: Λ is a datum; the minimal framework predicts w ≡ −1, **shared with ΛCDM** (F1). If DESI's
preference for evolving dark energy (DR2: 3.1σ with CMB; 2.8–4.2σ with CMB and supernovae, compilation-dependent;
no combination above 4.2σ) reaches a model-independent ≥5σ, the minimal framework and ΛCDM fall together and an
additional light field enters P0's content with a price (`03`: Cassini and Eöt-Wash).

## O3 — The outcome rule (corrected, R5)

P5 (non-contextuality) buys Born via Gleason/Busch. The revision-1 route via Hughston–Jozsa–Wootters is withdrawn:
HJW presupposes the measurement framework. The relevant literature: Gleason (1957); Busch (2003); Zurek's envariance
(2003–05); Deutsch (1999) and Wallace (2012) decision-theoretic derivations; Masanes–Galley–Müller (2019), who derive
POVMs, the Born rule and state update from unitary quantum mechanics plus operational closure, an associative
composition of outcomes and finite-dimensional state estimation — without assuming non-contextuality, a stronger
result than P5's Gleason route. Each trades the postulate for other assumptions; none derives probabilities from
causality alone. **GRUT III adds nothing here; P5 stays a postulate.** No computation proposed.

## O4 — The arrow (P4)

Located, not derived: the choice of the future horizon as the traced-out one. The physics is the cosmological initial
state. No computation proposed.

## O5 — The residual program: the only place a distinctive result could come from (new; owner's note, R11)

**The rule.** In any KMS state N_f = ½coth(βω/2)ρ_f (`07` T1), so GRUT III is standard wherever its principle
applies. A distinctive residual must be computed in a state where KMS fails, as Δρ_f (or ΔW(t,t′) for non-stationary
states, or Δκ_n) against the conventional baseline for the same observable, state and preparation. (The record's
G2-12 "Scope-Q" ruling, as the owner reports it, reached the same conclusion from the other direction: distinctive
locks need KMS-failing platforms. That ruling is not in this repository; the convergence is noted.)

**Baseline computed (`v11`, D11).** Quasi-de Sitter, conformal matter: the O(ε) KMS defect is fixed completely —
an effective-temperature shift T_eff ≈ (H/2π)(1 − 1.01ε) plus a 1%-of-ε frequency dependence. No free function.

**What the constructive branch's CP-05 settles for this program** (branch `codex/grut-causal-response-05`; review
relayed by the owner): one specified local scalar, eliminated, produces a covariant gravitational memory with derived
static range, slip and pole — a bounded physical answer, classified conventional scalar–tensor EFT. By `07` T6 the
kernel does not select that scalar; by `07` §4 the mechanism is an added field with a displayed price (`03`). Any
future construction's novelty claim must therefore be a *selection* statement forced by P0–P5, tested by connecting
independently measured quantities without refitting — not another realization of a passive kernel.

**Next computations, in order.**
1. The gravitational sector's own quasi-dS baseline: the massless (ν = 3/2 + O(ε)) worldline two-time function in
   power-law inflation, its stationary part's KMS defect, and the O(ε) correction to the zero-mode diffusion
   constant D = H³/4π² (expected: D(t) = H(t)³/4π² with the spectral tilt n_T = −2ε encoded in the two-time
   structure). This fixes what conventional theory says for the free graviton sector.
2. The one-loop conformal-matter correction to the tensor two-time function in quasi-dS (Fröb–Roura–Verdaguer did
   exact dS): the first *interacting* baseline.
3. Only then: is there any ΔW(t,t′) in the interacting gravitational infrared that conventional stochastic gravity
   leaves undetermined (the secular-growth problem)? If the literature's resummations disagree among themselves, that
   disagreement is the exact size of the opening; a GRUT statement would be a *restriction* on it forced by P0–P5.
   If they agree, the opening is closed and GRUT III is, finally and honestly, GR + QFT with a thermodynamic reading.

## O6 — a₀ (Milgrom's conjecture, not GRUT's; `06`)

Dead as a Jacobson-route argument (κ cancels). Alive only as Milgrom's vacuum-effect proposal with the
Narnhofer–Peter–Thirring/Deser–Levin temperature. Not pursued within GRUT III.

## O7 — Owed independent verifications of revision 2

1. Equation-level comparison of D2's worldline polynomial with Pérez-Nadal–Roura–Verdaguer (2009) and Cho–Hu (2014)
   (`08` records what could be read).
2. The massive-vector threshold polynomial for D5 (not computed; HLZ/GRW primary texts).
3. Independent reruns and adversarial review of `v7`–`v11`: the scheduled three-lens pass on v7–v10 **did not run**
   (usage limit); all five are author-level, with the author's own re-derivations recorded in `08`. First task of the next pass.
4. The record correction in D1 by re-running `calc/worldline_reduction.py` with the one-sided/symmetrised
   distinction made explicit.
5. The O(ε) defect of D11 by an independent method (e.g. adiabatic expansion of the conformal vacuum) and for the
   massless sector (O5 item 1).

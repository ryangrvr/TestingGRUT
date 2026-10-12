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

**The rule.** In any KMS state N_f = ½coth(βω/2)ρ_f (`07` T1), and P0's operator content is conventional, so ρ_f is the
conventional spectrum; together these make GRUT III standard at the two-point level wherever its principle applies (the
FDT alone fixes N given ρ; the second clause is what fixes ρ). A distinctive residual must be computed in a state where KMS fails, as Δρ_f (or ΔW(t,t′) for non-stationary
states, or Δκ_n) against the conventional baseline for the same observable, state and preparation. (The record's
G2-12 "Scope-Q" ruling, as the owner reports it, reached the same conclusion from the other direction: distinctive
locks need KMS-failing platforms. That ruling is not in this repository; the convergence is noted.)

**Baseline computed and corrected (`v11` rev. 3, D11).** Quasi-de Sitter, conformal matter: to first order the state
is *exactly* thermal at T_eff = (ä/ȧ)/2π = (H/2π)(1 − ε) for any a(t) with ε and η = ε̇/Hε first-order small (the revision-2 claim of a first-order frequency
dependence was a finite-ε artefact, caught by two reviewers); the first KMS violation is O(ε²) with the computed
shape D₂(ω). No free function. Today ε₀ = 0.47, η₀ = −2.1: the first-order form does not apply (D11; the ΛCDM kernel
differs from it by 0.1–19% over Δ = 0.25–1.5/H₀ in a grid-limited time-domain comparison that is neither a KMS test nor a
temperature scan, `11`), and only the exact kernel may be used for any late-universe statement.

**What the constructive branch's CP-05 settles for this program** (branch `codex/grut-causal-response-05`; review
relayed by the owner): one specified local scalar, eliminated, produces a covariant gravitational memory with derived
static range, slip and pole — a bounded physical answer, classified conventional scalar–tensor EFT. By `07` T6 the
kernel does not select that scalar; by `07` §4 the mechanism is an added field with a displayed price (`03`). Any
future construction's novelty claim must therefore be a *selection* statement forced by P0–P5, tested by connecting
independently measured quantities without refitting — not another realization of a passive kernel.

**A first-order question answered (`09`, `v12`).** The comoving detector temperature (ä/ȧ)/2π and the apparent-horizon
Kodama–Hayward temperature (H/2π)(1 − ε/2) differ at O(ε); P1 says only the former is a bath temperature (the apparent
horizon is not Killing for ε ≠ 0, so no KMS parameter exists for it), and the horizon first laws are identities of D3's
Einstein equations, not inputs. No heat-flow law; no O(ε) non-thermal window for the conformal sector.

**Literature verdict (`10`, fourth round).** Item 3: the scalar-sector resummations agree (stochastic = leading log, NLO
consistent) — closed; the graviton-sector dispute is live but is about which quantities are observables in perturbative
quantum gravity on dS, which P0–P5 do not address — not a GRUT opening. Item 1 is published except for the O(ε)
worldline KMS defect of the massless sector; its O(ε) zero-mode Fokker–Planck coefficients are Cho–Kim–Kitamoto (2015).

**The decisive cheap computation, attempted and restated (`v13` rev. 2; fifth and sixth rounds).** Question: does P1's
thermal reading at any horizon temperature fix the published O(ε) Fokker–Planck coefficients of the massless minimally
coupled zero mode in power-law FRW (Cho–Kim–Kitamoto 2015, eqs. 3.14–3.16)? Rev. 1 took the published coefficients, found
that the implied zero-mode temperature T_zm = (3−ε)A₀H⁴/2 (= 1/Vol(S⁴) at ε = 0) has O(ε) coefficient 2/3 − 2γ_E − 2 ln 2,
and concluded that no horizon temperature reproduces it (required c_X = 0.4685; detector 1, Kodama–Hayward 1/2, naive 0).
The reviewer showed the inference is not established: the stochastic noise is evaluated at the coarse-graining scale
k = σaH and σ is a choice. `v13` rev. 2 derives from the Bunch–Davies mode function that A₀(σ) = A₀(1)·σ^{−2ε/(1−ε)},
σ-independent in exact dS (Starobinsky) and not at O(ε): c₁(σ) = 1 − 2γ_E − 2 ln(2σ), c_X(σ) = 0.4685 + ½ ln σ. So c_X = 0 at
σ = 0.39, = 1/2 at σ = 1.07, = 1 at σ = 2.9, and c_X → −∞ as σ → 0 (a formal limit: the expansion needs ε|ln σ| ≪ 1). The γ_E and ln 2 of
rev. 1 are γ_E + ln(2σ) at σ = 1: a matching-scale artifact. **The sharp-window O(ε) zero-mode coefficient is scheme-dependent and cannot test horizon temperatures; rev. 1's
inference is withdrawn** (its arithmetic at σ = 1 stands). Two corrections from the sixth-round audit (`11`): (i) a σ-free
invariant *does* exist — the zero-frequency noise strength D(t′) = lim_{t→∞}∂_{t′}W(t,t′), with the closed form
D/H³ = (√π/2π³)(1−ε)³Γ(ν)Γ(5/2−ν)cos(πε/(1−ε)), whose O(ε) coefficient is exactly −(1 + 2 ln 2), the number `v14` E finds from
the two-time function (`v13` 13.5 checks the coefficient symbolically). So "any dimensionless comparison shifts by a
logarithm" was too strong: a diffusion coefficient is invariant; a *temperature* read off the zero mode is not. (ii) The
zero-mode-temperature route fails for a prior reason: Cho–Kim–Kitamoto's own late-time solution for ε ≠ 0 (their eq. 3.27) is
exp(−V/T_zm) times an extra O(ε) Gaussian factor, so no stationary weight exp(−V/T) exists at O(ε); and their drift
1/((3−ε)H) is a leading-IR approximation (the exact k = 0 attractor gives 1/((3+ε)H)), so the "published conventional
answer" is itself approximate at this order. The σ-independent object for the *temperature* question is
frequency-resolved: the comoving detector's Wigner spectrum at ω ≠ 0 and its KMS ratio at O(ε), the quantity for which the
conformal sector gave exact thermality at (ä/ȧ)/2π (`v11`). That computation is `v14`, below. **What stands:** O5 items 1
and 3 close on other grounds — P1 assigns no temperature when the modular flow is not geometric, and the σ-independent
correlators are fixed by conventional QFT.

**The σ-free computation, done (`v14`, sixth round; D13).** The massless minimally coupled sector's comoving two-point
function in power-law FRW has a closed form (the de Sitter function of index 3/2 + ε/(1−ε) times a ratio of scale factors,
checked against a direct mode integral). At first order in ε it is (1 + 3ε/2)·W_mmc^{dS}(Δ; ä/ȧ) + ε(H̄²/8π²)·S₂(cosh²(H̄Δ/2))
+ ω = 0 terms, where S₂ is an explicit dilogarithm. Because any function of cosh²(H̄Δ/2) is KMS at 2π/H̄ by itself, the
first-order KMS defect of the ω ≠ 0 part at T = (ä/ȧ)/2π vanishes (checked on the exact closed form): **the fixed-mean-time
(Wigner) spectrum of the massless sector is thermal at (ä/ȧ)/2π to O(ε), exactly like the conformal sector's, and the
Kodama–Hayward rate is excluded at O(ε)** — no coarse-graining scale enters. (Convention, from the audit: this is the Wigner
spectrum, not an instantaneous one-sided detector rate, which differs at O(ε) by a frequency-dependent Ḣ term and, for this
sector, is non-thermal already in exact dS unless derivative-coupled; D13.) The first-order effect is a shape change of the spectral density, ρ = ρ_mmc(ω; H_eff)[1 + ε r(ω)] with
r(0) = 2 − 2 ln 2 and r(∞) = 0, and the σ-free zero-mode diffusion coefficient is D = (H̄³/4π²)[1 − (1 + 2 ln 2)ε], which the
stochastic formula reproduces at σ* = e^{1−γ_E} = 1.53 (a calibration). **O5 item 1 is now computed, not argued** — for
constant ε and in the Wigner convention; the general slow-roll extension is open (D13 Scope). It is conventional QFT and
establishes nothing about P1 beyond consistency.

**What remains, in order (after `v13`).**
1. (closed) Items 1 and 3 — the massless-sector baseline (`v14`: thermal at (ä/ȧ)/2π at first order; D13) and the
   secular-IR opening (`10`).
2. The one-loop conformal-matter correction to the tensor two-time function in quasi-dS (Fröb–Roura–Verdaguer 2012 did
   exact dS): the last conventional baseline. Low priority: by the precedent of `v11` (conformal sector, O(ε²)) and `v14`
   (massless sector, O(ε)), conventional theory fixes it and P1 adds nothing.
3. The graviton infrared (which observables carry secular logs) is a question in perturbative quantum gravity on dS;
   none of P0–P5 addresses it, and GRUT III does not claim it.
4. Independent verification of `v8`, the general-a(t) and ΛCDM parts of `v11`, `v12`, `v13` (O7).

**The honest closing line now applies.** GRUT III is GR + QFT with a thermodynamic reading. In every regime computed in
this branch — exact dS (KMS by construction), quasi-dS conformal sector through O(ε²) (`v11`), quasi-dS massless sector
at O(ε) (`v14`, σ-free), the horizon first laws (`v12`) — the answer is fixed by conventional theory and P0–P5 add no selection.
The program's genuine fork is therefore not another computation inside P0–P5 but a decision: either (a) GRUT III is
accepted as the thermodynamic *reading* of GR + QFT, a consistent framework with no distinctive prediction (`03`), or
(b) a sixth postulate is proposed that breaks KMS-by-construction in a controlled, parameter-free way — the only route to
a distinctive prediction — and is tested with the Δρ instrument of `07` against the baselines fixed here, priced against
the fifth-force bounds of `03`. Option (b) is the open road; nothing in this branch yet supplies such a postulate. A sixth
postulate of that kind is exactly a Stage-3 𝒦 card: the frozen charter on `grut2-stage3` (c1275a8) already supplies the
pricing, the relocation tests, the Scope-Q non-KMS requirement and the price list such a card must meet. If option (b) is
taken, it is written as Card 1 under that charter — after the owner review its hard stop (§25, G2-11) requires — not as new
governance here. That merges the two tracks instead of rebuilding one.

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

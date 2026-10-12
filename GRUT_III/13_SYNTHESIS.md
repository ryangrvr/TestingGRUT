# 13 — Synthesis of all the guidance: what is right, what cannot be right, and the one open door

**Status (2026-10-12).** Written at the owner's request to synthesise the guidance received so far, on the assumption that
some of it is wrong. The guidance: the GRUT III baseline and its Card-1 constraint sheet (`12`); the reviewer's fork
recommendation; ChatGPT's self-assessment of the program; and the codex graph line (fixed-bond entropy theorem, SCOUT-RG01,
SCOUT-RG02, CP-12 and its three audits, read from their branches). Method: every claim the synthesis rests on was checked
against primary sources by a six-group literature workflow and a codex-mathematics audit (§8 records verdicts and read
levels). Two new computations were added: `v15` (the Sorkin–Johnston state in a de Sitter slab) and `v16` (the
entropy–action budget). Where the checks contradicted my own earlier guidance, it is corrected here and in `09` and `12`.

## 0. The answer in one paragraph

Nothing in either line has found new physics, and the reasons are now understood rather than felt. GRUT III reads the
temperature off the state's own modular flow, so it can add nothing to the state that quantum field theory already fixes.
The graph line's thermal route is excluded by counting: Gibbs weights at a fixed positive temperature cannot make geometry
dominate a configuration space whose non-geometric part is super-exponentially larger, as N → ∞. `v16` argues,
heuristically, that a generic local action's phase cannot either. Its zero-temperature ground states escape that argument
but failed for model-specific reasons (built-in valence, fragmentation with N, no three-dimensional or Lorentzian
structure). The graph line's target, bounded-degree spatial graphs with a Lieb–Robinson cone in an external time, would
have a microscopic preferred frame. Unless a protecting mechanism is supplied, generic loop effects then make species'
limiting speeds differ by 0.1–10% (Collins et al. 2004), against one-sided astrophysical bounds of ~10⁻²⁰. What survives is
two narrow doors where local covariance leaves physics a genuine choice: the quantum *state*, and the finitely many
*renormalisation constants* of the stress tensor, among them Λ and G. A distinctive GRUT prediction has to come through one of
them, fixed by a principle rather than a temperature or a fit. For the state, every known rule fails or needs a scale the
principles do not fix. The candidate free of continuous parameters, the Sorkin–Johnston state, is computed here (`v15`): it
violates GRUT's own P2 in a sharp region and contradicts the CMB over the inflationary slab. The unresolved freedom sits in
the infrared of massless fields in an accelerating universe. That is where the state is least determined, where GRUT's
horizon principle could act, and where observation is weakest.

## 1. What each source of guidance said, and its verdict

- **GRUT III baseline (`12`)**: a verified thermodynamic reading of GR + QFT, no distinctive prediction. *Stands.* Two of its
  lines were overstated and are corrected below (de Sitter KMS needs Hadamard; P2's state content is the Hadamard condition).
- **The reviewer**: bank the baseline and write the sixth postulate as Card 1 under the frozen charter; the options are not
  exclusive. *Right.* This synthesis turns "a sixth postulate that breaks KMS in a controlled way" into two precise targets
  (§5: a state-selection rule, or a principle fixing a renormalisation constant) and adds three admissibility lines (`12`
  constraints 8–10).
- **ChatGPT's self-assessment**: much mathematics, no new law, progress over-reported; demand originality audits, constrained
  construction and an external benchmark. *Right in substance.* Its proposed replacement principle ("derive the configuration
  space from constraints") points the right way but is too vague to act on; §3 gives the precise requirement. The external
  benchmark it asked for exists and is decisive: Lorentz-invariance data (§3).
- **The codex graph line**: correct model-specific theorems; no spacetime. *Its mathematics stands (§8); its premise cannot be
  right* (§3).
- **My own earlier claim that "no state-selection rule has been proposed"** (`12` §3): *wrong.* Several exist (Sorkin–Johnston,
  states of low energy, horizon states, adiabatic vacua). Each buys uniqueness with a cost; §5 prices them. The candidate space
  also missed a fourth route — fixing a renormalisation constant such as Λ by principle — now added to `12` §3.

## 2. The diagnosis: two lines, two different failures, one lesson

- **GRUT III, with geometry given.** P1 identifies the bath temperature with the KMS parameter of the state's own modular
  flow wherever that flow is geometric. That is a definition, not an extra input, and quantum field theory already fixes
  every correlator of that state. So the reading adds nothing: KMS by construction, no distinctive prediction (`12`;
  `v11`, `v14`).
- **The graph line, with geometry sought.** Its thermal route is excluded by counting. A Gibbs state at a fixed positive
  temperature, over a configuration space whose non-geometric part outnumbers the geometric part super-exponentially,
  cannot select geometry with extensive energies as N → ∞. Its zero-temperature route is not covered by that theorem;
  it failed for model-specific reasons (an exact T = 0 selection result that rests on built-in valence, quantum ground
  states that fragment as N grows, no three-dimensional or Lorentzian structure).

The shared lesson is narrower than a slogan. Thermality cannot be the *source* of anything: once a geometry and a state are
given it is a property of the state, and before there is a geometry there is no bath — the universe has no outside — so
"a fixed positive temperature over pregeometries" is not a physical premise. In Lorentzian approaches positive weights are
not what quantum gravity sums; Euclidean and Wick-rotated approaches do use them, at the price of a continuation back. Two
corollaries organise everything below:
1. Geometry has to come from **kinematics plus amplitudes**, not from thermal selection: a configuration space small enough
   (or an action large enough) that the weights can single out manifold-like structure (§3).
2. A distinctive prediction has to come through one of the two choices local covariance leaves open: **the state**
   (Fewster–Verch 2012 forbid any covariant state assignment), or **the renormalisation constants** of the stress tensor
   (Λ, G and the curvature-squared coefficients, fixed by locally covariant renormalisation only up to these finitely many
   numbers; Hollands–Wald). GRUT's horizon relationality (P1) could act on the first; a horizon principle that fixed Λ would
   act on the second (§5).

## 3. The graph line: what is right, what cannot be right

**Right, and worth keeping.** The codex theorems are correct as stated (independently re-derived, §8), with their scope
limits. Conserving the number of connections removes the dense-edge failure. But at a fixed positive temperature a family
of lattice templates (at most N!·e^{O(N)} labelled members) has probability at most exp[−(ν−1)N ln N + O(N)] among graphs
with Q = νN connections. This holds even with the full non-commuting gauge–fermion–rewiring Hamiltonian (relative-entropy
bound). The CP-12 back-reaction (a charged connection is removed at order T⁴ instead of T²), the Gauss-law distance
diagnostic, the O(T²/N) operational window and the coherent no-assembly bound are sound model-specific results. Two
qualifications from the audit: the CP-12 thermal exclusion is an O(1/ln N) upper bound that is vacuous below N ≈ 10⁷–10⁸ at
the frozen couplings, and the operational window fails for coherent superpositions of graphs.

**Rediscovered, not discovered.** The counting obstruction is decades old in other notation. Chen & Plotkin (Phys. Rev. D 87,
084011, 2013) found numerically that the transition temperature of fixed-edge graph models goes to zero as N grows.
Kleitman & Rothschild (1975) is the same kind of obstruction for causal orders. Quantum graphity (Konopka, Markopoulou &
Severini 2008; Hamma et al. 2010) met the locality and equilibrium problems the RG02 and CP-12 audits found. Trugenberger and
Kelly–Trugenberger–Biancalana obtained geometric phases only by restricting the configuration space by hand (regular
bipartite graphs) and rescaling the coupling with N — the two moves the codex charter forbids. Euclidean dynamical
triangulations are the instructive contrast: at fixed topology their counts are exponentially bounded, so they escape the
counting obstruction, but they met a different one — their phases are crumpled or branched-polymer rather than extended 4D
geometry — which causal triangulations later addressed with a foliation.

**The precise reason (`v16`).** Selection needs an *entropy–action budget*. A weight can make a geometric class dominate
only if what it can vary across the competing class is at least that class's excess log-count.
- Boltzmann weights with extensive energies (|E| ≤ CN), at a fixed positive temperature: the codex theorem — the excess
  (α−1)N ln N wins as N → ∞.
- Lorentzian weights e^{iλS} with a local action on sparse graphs: the class sum is a product of O(N) roughly independent
  local factors, so generically it is suppressed by exp(−N g(λ)) with g = O(1). At weak coupling g ≈ λ²Var(S)/2N (`v16`
  checks this Gaussian regime at λ = 0.05, with Var S ∝ N). At strong coupling g is smaller and can vanish: for the degree
  action S is always even, so at λ = π there is no suppression at all. Either way it is exp(−O(N)) or none, and that loses
  to N ln N. This is a heuristic, not a theorem: exact cancellations at tuned couplings are not excluded, and the phase also
  acts on the geometric class. Interference is therefore not expected to rescue local actions.
- The two ways out that are known: a configuration space whose counts grow only exponentially (triangulations of fixed
  topology: proven in 2D, for locally constructible 3D triangulations and for 3D causal triangulations; open in general 3D
  and in 4D), or a nonlocal action with ~N² terms against N² entropy (causal-set theory: the dominant layered orders are
  suppressed by the link or BDG action — Loomis–Carlip 2018 for two-level orders; Mathur–Singh–Surya 2021 and
  Carlip–Carlip–Surya 2023–24 for Kleitman–Rothschild/layered orders — at leading order and in periodic parameter bands).
  No paper yet compares these sums with the manifold-like contribution, and the full two-level class is not exponentially
  small, so this is necessary and promising, not solved.
- One sharpening of the codex result (graph-model check, §8): if vertex relabellings are treated as gauge (unlabelled
  counting), the competing class is super-exponential for *every* ν > 0 — at least exp(νN ln N − 2νN ln ln N + O(N)) — so
  CP-12's "ν = 1 is open" closes in that reading. In labelled counting it stays open.

**The target owes an explanation that no known mechanism supplies.** CP-12's stated goal is a bounded-degree, spatially
organised graph on which matter propagates inside a Lieb–Robinson cone in external time. Two facts and a benchmark:
1. *An external clock picks a frame.* A spatial graph evolved by a Hamiltonian in an external time has a rest frame by
   construction, and its Lieb–Robinson cone is a model-dependent upper-bound cone in that frame. For discreteness meant to be
   Lorentz invariant, Bombelli, Henson & Sorkin (Mod. Phys. Lett. A 24, 2579, 2009) prove there is no Lorentz-equivariant
   measurable rule assigning a preferred direction to a Poisson sprinkling of Minkowski space; that no finite-valency graph
   can be attached equivariantly is argued informally in the paper. The theorem does not cover abstract, unembedded graphs —
   for those the frame comes from the clock.
2. *A frame at the Planck scale does not stay there, unless something protects it.* Collins, Perez, Sudarsky, Urrutia &
   Vucetich (PRL 93, 191301, 2004): in an effective theory with a Lorentz-violating ultraviolet regulator, ordinary loops
   generate unsuppressed dimension-4 terms that shift each species' limiting speed by c → c(1 + ξ/4), ξ = (g²/6π²)[1 + …] —
   0.1% to 10% for Standard Model couplings, independent of the cutoff scale. In an emergent model species' speeds need not
   agree even at tree level, so the problem could be worse, or absent if all species share one origin. Proposed escapes —
   low-scale supersymmetry (Groot Nibbelink–Pospelov 2005; not observed), strong-coupling power-law running
   (Bednik–Pujolàs–Sibiryakov 2013), gravitational confinement of the violation (Pospelov–Shang 2012) — are proofs of
   principle without a Standard Model embedding; logarithmic infrared running (Chadha–Nielsen 1983) is far too slow.
3. *The benchmark.* Limiting-speed differences between species are bounded at ~10⁻²⁰–10⁻²¹ (electron–photon,
   proton–photon; one-sided astrophysical bounds; two-sided laboratory bounds reach ~10⁻¹⁵), photon birefringence at ~10⁻³⁵,
   and linear energy-dependent photon speeds above ~2 E_Planck (Vasileiou et al. 2013, conservative; Fermi GRB 090510:
   > 1.2 E_Pl) (Kostelecký–Russell tables, 2026 edition). A bounded-degree pregeometry with an external clock owes an
   explanation of that equality. This is the external, unfitted benchmark the codex line was asked to find.

**External time.** Every codex model (CP-11, CP-12, RG01, RG02) evolves in an external parameter with a Hamiltonian bounded
below, and none has a causal order as a primitive (§8). A distinguished time is not by itself excluded: unimodular gravity
and shape dynamics are equivalent to general relativity with one, and causal dynamical triangulations use a foliation. But
a cone that is a Lieb–Robinson bound in that time is the same Lorentz-invariance debt as items 1–3, seen from the dynamics.

**Verdict on the line.** Its obstructions are real and now understood. Its premise — undirected adjacency, Gibbs weights,
external time, bounded degree — has each of its four ingredients unexplained, and two of them (Gibbs weights at a fixed
temperature; bounded degree with an external clock) carry known obstructions. Continuing to add mechanisms inside that
premise will keep producing correct theorems about models that do not yet describe our universe. If the line continues, it
should change kinematics to a causal order (directed, acyclic) plus a counting measure, with Lorentzian amplitudes and an
action whose budget matches the entropy. That is causal-set theory's choice, made for exactly these reasons (§7 states its
open problems); GRUT should build on it rather than re-derive its obstacles.

## 4. Corrections to the GRUT III guidance (mine)

- **De Sitter KMS needs Hadamard.** `12` said the de Sitter vacuum is KMS "by de Sitter invariance and positivity". De Sitter
  invariance alone allows the Mottola–Allen α-vacua, which are not thermal for geodesic observers. Bunch–Davies is the unique
  de Sitter-invariant pure Gaussian *Hadamard* state for m² > 0, and for the massless minimally coupled field no de
  Sitter-invariant state exists at all (Allen 1985; Hollands–Wald 2015). Constraint 1 of `12` is corrected.
- **P2's state content is the Hadamard condition, as a scaling-limit statement.** Jacobson's 1995 "local equilibrium" has a
  geometric part (an instantaneously stationary local Rindler horizon) and a state part ("the usual short-distance vacuum
  fluctuations"). Only the state part is the Hadamard condition. Exact KMS at κ/2π follows for quasifree, Killing-invariant states that
  are Hadamard near a bifurcate Killing horizon with a wedge reflection (Kay–Wald 1991), and such states need not exist (Kerr,
  Schwarzschild–de Sitter); at other horizons it holds in a near-horizon scaling limit (Kurpicz–Pinamonti–Verch 2021;
  Moretti–Pinamonti 2012). The identification has prior art: Chirco, Haggard, Riello & Rovelli 2014. Added to `12` as
  constraint 8.
- **The apparent-horizon temperature is more than a label at the horizon.** `09` called T_KH a geometric label. Kurpicz,
  Pinamonti & Verch prove that Hadamard states become KMS with respect to the Kodama flow at the Hayward surface gravity, in
  a near-horizon scaling limit, for spherically symmetric outer trapping horizons. T_KH and T_det = (ä/ȧ)/2π are then two
  temperatures of two different observables, both conventional. Whether the result covers the cosmological apparent horizon
  is open. `09` is amended.
- **The "no natural state" theorem is Fewster–Verch 2012** (Ann. Henri Poincaré 13, 1613, Thm 6.13), not
  Brunetti–Fredenhagen–Verch 2003, which only sketched it. It forbids one state assignment compatible with *all* spacetime
  embeddings for a dynamically local theory. It does not forbid distinguished states in spacetimes with extra structure. That
  is exactly why a selection principle has to use global or structural data.

## 5. The two places a principle could act — and the test of one existing candidate (`v15`)

Every distinctive claim has to change the operator content, the identification, the state, or a constant that local
covariance leaves free (`12` §3, now with four routes). Changing the operator content is conventional physics with a price.
Changing the identification abandons GRUT III. That leaves two slots.

**The second door: renormalisation constants.** Locally covariant renormalisation fixes the stress tensor only up to
finitely many local curvature terms, whose coefficients include Λ and G (Hollands–Wald 2001–2005). A principle that fixed
one of those numbers — Λ from a horizon principle, for instance — would be distinctive without touching the state. No such
principle exists on this branch (`06`'s a₀ ∼ cH₀ is Milgrom's coincidence, not a derivation). It is listed as N2 in §9.

**The first door: the state.** Here the field has a precise, known gap. Fewster & Verch (2012, Thm 6.13) prove that no state
assignment can be compatible with all spacetime embeddings for a dynamically local theory. Distinguished states exist only
where extra structure is supplied. In de Sitter the Bunch–Davies state is selected by symmetry plus Hadamard. In FRW every
known rule buys uniqueness with an extra input: an initial time (adiabatic vacua), a smearing function (Olbermann's states
of low energy), a smoothing function (Brum–Fredenhagen 2014), a region (Sorkin–Johnston), or a past horizon with
past-asymptotic de Sitter (Dappiaggi–Moretti–Pinamonti). `v13` showed what a hidden scale does: it can move a first-order
coefficient anywhere. A principle that selected the state from structure the universe actually has would be a genuine
addition. This is where P1 (the bath is the observer's causal complement) naturally speaks, because a rule that local
covariance cannot supply has to use global data.

**One candidate, computed.** The Sorkin–Johnston (SJ) prescription (Johnston 2009; Afshordi, Aslanbeigi & Sorkin 2012) needs
no Killing field, no temperature and no test function; its one input is the region. It keeps the positive part of the
Pauli–Jordan operator on the region's L² space. `v15` reduces it exactly, mode by mode, in a de Sitter slab ending at a_eH:
two integrals per mode, with the SJ condition ∫v_k² a⁴dη = 0. The reduction agrees with the positive part of the discretised
operator to 4e-12 for the conformal field at k = a_eH, and to six digits for the massless field deep in the infrared. An
independent critic re-derived every result analytically. It finds:
- **Ultraviolet:** |β_k|² → (a_eH)²/(16k²) for both the conformal and the massless field, from the boundary term at the slab's
  end. That is a power law, so the state is not Hadamard (a direct computation; Fewster–Verch 2012 proved the analogous result
  for massive fields on compact ultrastatic slabs). The |β|² part of its energy density grows like Λ² with the cutoff. It
  violates P2.
- **Infrared:** modes outside the horizon at the region's end are suppressed. To leading order for a_iH ≪ k ≪ a_eH,
  R_k = |v_k|²/|u_k|² ≈ √(π/2)(k/a_eH)^{1/2} for the conformal field and ≈ √(π/2)(k/a_eH)^{3/2} for the massless field. The
  prefactors come from ∫sin²x/x² = π/2 and ∫j₁² = π/6, and power counting gives (k/a_eH)^ν for a light field. If the region is
  the inflationary slab, every observable mode has its power multiplied by ~10⁻³³–10⁻³⁹ and its tilt shifted by +ν ≈ +1.5,
  against Planck's n_s − 1 = −0.035 ± 0.004. Excluded in amplitude and tilt, for any light inflaton.
- **Region dependence:** moving the region's future boundary from a_eH = 1 to 100 changes the two-point function of a
  k = 10⁻³ mode on the common earlier region by a factor of about 10³. This is the known non-local covariance of SJ
  (Fewster–Verch 2012 §5; Afshordi–Aslanbeigi–Sorkin 2012 §6). A rule that is causal in GRUT's sense must not depend on a
  future boundary.
- **Not computed:** a region ending today. That needs the reduction through reheating and the radiation and matter eras; the
  UV tail would then sit at (a₀H₀)²/(16k²). Whether it is compatible with the large-angle CMB is open.

**Where this sits in the literature.** Where it is defined, the limiting SJ state in de Sitter is an α-vacuum
(Aslanbeigi–Buck 2013). On the global patch in 4D it is Bunch–Davies for m ≥ 3/(2ℓ) and at the isolated mass mℓ = √5/2. For
other lighter fields it is a non-Hadamard α-vacuum or ill-defined (including the conformally coupled massless field); on the
Poincaré patch it is undefined below that mass. The finite-slab results of `v15` for the light fields are the finite-region
face of the same facts; the explicit UV tail and IR exponents are direct computations, not quoted. Afshordi, Aslanbeigi &
Sorkin (2012) found that in a radiation-era FRW the massless SJ state is the conformal vacuum, and they stress the region
dependence.

**The other unique rule.** The horizon state of Dappiaggi, Moretti & Pinamonti (2009) comes closest to GRUT's P1. It is built
from data on the past cosmological horizon, is unique and invariant under the horizon's symmetry group, equals Bunch–Davies
in exact de Sitter, and is proven Hadamard for spatially flat FRW that is asymptotically de Sitter in the past. **But the
proof needs Re ν < 3/2**, with ν = (9/4 − m²/H² − 12ξ)^{1/2} fixed by the asymptotic past Hubble rate (or Re ν < 1/2 under
weaker asymptotics). It therefore excludes the massless minimally coupled field (ν = 3/2) even in backgrounds that are
asymptotically de Sitter in the past. A slow-roll background with ε bounded away from zero in the past is not asymptotically
de Sitter at all; there the field's mode index is 3/2 + ε/(1−ε) (the `v12` lemma). So the GRUT-native rule that works is
not established for the fields that seed cosmic structure.

**Summary of the state door.** The best-studied region-based rule, SJ, fails GRUT III's own P2 in a sharp region and, over
the inflationary slab, contradicts the CMB; smoothing it restores Hadamard and reintroduces a free function. The horizon
rule avoids the region choice but is not proven at ν = 3/2. The target for a Card 1 through this door is a rule that is
Hadamard, causal (no future boundary), free of arbitrary functions, and that reproduces Bunch–Davies for the modes we observe.

## 6. Most routes end at a free scale, and the unresolved ones in the infrared

Every attempt in this program to get beyond the conventional answer has stopped at a choice the principles do not fix.
The unresolved ones are infrared choices for massless fields in an accelerating spacetime.

| Route | Where it stopped |
|---|---|
| `v13` zero-mode "temperature" | the coarse-graining scale σ moves the O(ε) coefficient anywhere |
| `v14` massless sector | the IR constant and the secular term, which the Wigner convention puts at ω = 0 |
| `10` secular-IR survey | the graviton dispute is over which IR quantities are observable |
| Sorkin–Johnston (`v15`) | UV: non-Hadamard (power-law tail). IR: the region's future boundary, through k/a_eH, decides the superhorizon state |
| low-energy and smoothed states | a temporal smearing or smoothing function (its effect is expected, not shown, to concentrate in modes with k/a below the inverse support time) |
| horizon state (DMP) | the Hadamard theorem needs Re ν < 3/2, the IR-integrability condition; the massless minimally coupled field (ν = 3/2) is excluded |

For the massless field this is no accident. For ν ≥ 3/2 a Bunch–Davies-like two-point function has the Γ(−δ) pole that `v14`
expanded around, and ∫k²|u_k|²dk diverges at small k — logarithmically at ν = 3/2, as a power beyond (Ford–Parker 1977;
Allen 1985). The state of those modes is fixed by when and how the acceleration began, not by local physics. So the infrared
is where conventional theory genuinely leaves the state open. It is where a horizon-relational principle like P1 could act.
And it is the one place observation cannot reach well: modes that left the horizon near the start of inflation are today at
or beyond the largest observable scales, where cosmic variance dominates. Planck finds every non-Bunch–Davies bispectrum
template consistent with zero, and any oscillatory modulation of the power spectrum below about 3%.

**Consequence.** A GRUT prediction through the state, if one exists, lives in the lowest multipoles of the CMB and the
largest scales of structure. It must be made with its scale fixed by the principle before looking, not fitted to the ~2σ
low-ℓ deficit. With the principles as they stand, the scale is not fixed, so there is no prediction.

## 7. What a theory built on these lessons would look like

Three layers, each with the benchmark it must meet. None is complete; this is the map, not the territory.

1. **Microstructure: causal order plus counting, not adjacency.** In the continuum, a causal order fixes a Lorentzian
   metric up to a conformal factor in d > 2 (Hawking–King–McCarthy 1976; Malament 1977), and volume supplies the rest.
   That a causal set determines its approximating spacetime is conjectured (the Hauptvermutung), not proven, and "number ~
   volume" holds only statistically. Sprinklings are statistically Lorentz invariant, so Lorentz invariance is built into the
   kinematics; whether the dynamics selects sprinkling-like orders is the open entropy and continuum problem, and no
   first-principles quantum dynamics exists. The price of the kinematics is nonlocality (in a sprinkling of all of
   Minkowski space, infinitely many links per element). *Benchmark:* manifold-like orders must dominate. Today the balanced
   layered classes are suppressed at leading order in periodic parameter bands; the result depends on the subset summed
   ("KR + dust" is not suppressed), and no comparison with manifold-like orders exists.
2. **Low-energy dynamics: general relativity from P2.** Jacobson 1995 and the entanglement-equilibrium version (Jacobson 2016)
   give the Einstein equation from local horizon thermodynamics, given a Hadamard-type state and a finite universal entropy
   density — cleanly for conformal fields at first order; for general fields Jacobson 2016 rests on a conjecture about the
   modular Hamiltonian, and finite variations are open. *Benchmark:* GR's tests, which GRUT III passes because it is GR.
3. **The state: Hadamard (P2), selected horizon-relationally (P1) — a programme, not a result.** Where a past cosmological
   horizon exists and the theorem applies (massive or non-minimally coupled fields, spacetimes asymptotically de Sitter in
   the past), the DMP state is distinguished and equals Bunch–Davies in de Sitter. *If* P1 can be shown to select it — not yet
   done — GRUT would recover the standard initial state from an assumption of past-asymptotic de Sitter rather than
   postulating Bunch–Davies directly. That would be explanatory, not predictive, and it would not yet cover the massless
   minimally coupled fields that seed structure (N1). *Benchmark:* Planck's non-Bunch–Davies and oscillation bounds, and n_s.

What GRUT adds, honestly stated: a candidate selection rule for layer 3 (unproven), a possible principle for a
renormalisation constant (none proposed), and an interpretation of layer 2. It adds nothing yet to layer 1, where
causal-set theory has the right kinematics and GRUT has no dynamics of its own.

## 8. Verification record

Seven verification agents (2026-10-12) checked every claim above against primary sources or by independent computation.
Read levels are given because this project has misattributed sources before. Full structured output is kept in the session
scratchpad; the verdicts that change or qualify the text are listed here.

**Codex mathematics (one agent, independent re-derivations and reruns of copies).**
- Fixed-bond Gibbs class bound, and its (α−1)N ln N consequence: *correct as written.* It needs an extensive operator-norm
  remainder, so it does not apply directly to unbounded rotors. It describes a counterfactual swap model, not frozen CP-11.
- CP-12 relative-entropy theorem, P_β = O(1/ln N) for ν > 1: *correct, but an upper bound set by the method, not a sharp
  rate.* At the frozen couplings it is vacuous (≥ 1) until N ≈ 10⁷–10⁸ at β = 1. It gives no exclusion for ν = 1.
- RG02 thermal single-ring probability ~ C_β N^{−1/2}: *consistent.* At β = 0 it is exactly the uniform perfect-matching
  fraction √π/(2√N); C_β = 0.886, 1.70, 3.63, 7.23 at β = 0, 1, 5, 20. It covers only the conditional ensemble.
- Uniform single-cycle fraction π e² (2/e²)^N: *correct.*
- Coherent-assembly limit 1 − exp[−4ν(1 − cos ΓT)] and the no-assembly bound O((1 + Γ²T²)/N): *correct.* They hold for the
  neutral comparator only, not for full CP-12 evolution.
- Operational window |R_full − R_fixed| ≤ 4Γ²T²Q(M−Q)/M²: *correct* for graph-diagonal sources and detectors and for
  graph-dephased states. It fails for coherent graph superpositions, which pick up O(T) terms.
- Structural survey: every model has an external time and a Hamiltonian bounded below. None has a causal order as a primitive.
  Every "cone" is a propagation bound in external time on supplied or prepared adjacency.
- Overstatements found: the committed `cp12_dynamic_distance_audit/check.py` crashes (only a repaired copy runs); the CP-12
  distance relation is for the bare electric energy (the full diagonal has slope (ε_E − a_N)/2); the RG02 detector hierarchy
  is an N ≤ 4 observation; the CP-12 Lieb–Robinson cone cites theorems whose constants were not computed. The largest CP-12
  eigensolver runs (up to 2.9 million states) were not rerun.

**Literature (six agents).**

| Claim used here | Verdict | Read level |
|---|---|---|
| Kleitman–Rothschild: almost all orders are 3-layer, log₂ count n²/4 + 3n/2 + O(log n) | confirmed | abstract + restatement |
| Loomis–Carlip suppress Kleitman–Rothschild orders | **refuted** — they suppress two-level orders; KR came later (Mathur–Singh–Surya 2021; Carlip–Carlip–Surya 2023–24) | full text |
| Suppression scales like the entropy (N² against N²) | confirmed for the balanced subclass; the full two-level sum is not exponentially small | full text + computed |
| Bombelli–Henson–Sorkin: no Lorentz-equivariant direction for sprinklings; corollary for finite-valency graphs | confirmed (corollary argued informally in the paper) | full text |
| Collins et al. 2004: Planck-scale frame gives O(0.1–10%) species speed differences | confirmed, with its assumptions (explicit Lorentz-violating regulator, no custodial symmetry) | full text |
| Species speed-difference bounds ~10⁻²⁰–10⁻²¹ (one-sided astrophysical), ~10⁻¹⁵ (two-sided lab) | confirmed, Kostelecký–Russell 2026 edition | full text |
| Fermi GRB 090510: E_QG,1 > 1.2 E_Pl; Vasileiou et al.: robust ≳ 2 E_Pl | confirmed | full text |
| "No natural state" | **re-attributed** to Fewster–Verch 2012, Thm 6.13 | full text |
| Sorkin–Johnston: not Hadamard on ultrastatic slabs (a.e. width), not locally covariant | confirmed | full text |
| SJ in de Sitter is an α-vacuum; Bunch–Davies only for m ≥ 3/(2ℓ) on the global patch | confirmed (Aslanbeigi–Buck 2013) | full text |
| DMP horizon state: unique, Bunch–Davies in dS, Hadamard for Re ν < 3/2 | confirmed; excludes the massless minimally coupled field | full text |
| Hadamard ⇔ microlocal spectrum condition | confirmed (Radzikowski 1996) | abstract |
| P2's state content is Hadamard | **partly correct**: as the formal statement of Jacobson's assumption, with KMS exact only at Killing horizons and otherwise in a scaling limit | full text |
| Bunch–Davies is the unique dS-invariant pure Gaussian Hadamard state (m² > 0); α-vacua not KMS | confirmed | full text |
| Planck: non-Bunch–Davies bispectrum templates consistent with zero; oscillation modulation ≲ 3% | confirmed | full text |
| Chen–Plotkin: T_c → 0 with N in fixed-edge graph models | confirmed (numerical, heuristic counting) | full text |
| Trugenberger / Kelly–Trugenberger–Biancalana: geometry by restricting to regular graphs and rescaling the coupling with N | confirmed | full text |
| Exponential bound on triangulations: 2D proven; 3D locally constructible and 3D causal proven; general 3D and 4D open | confirmed | full text |
| CDT 4D de Sitter phase | confirmed; the 2008 paper is Ambjørn, Görlich, Jurkiewicz & Loll | full text |
| Unlabelled graphs with αN edges ≥ exp((α−1)N ln N + O(N)); with all edges on few vertices ≥ exp(αN ln N − 2αN ln ln N) for every α > 0 | confirmed (computed) | computed |

Where a read level is "abstract", the claim rests on the abstract plus a full-text restatement in a later source.

**Adversarial review of this document (two critics, 2026-10-12).** A physics-and-framing critic and a numerics critic tried
to refute the draft. Neither found a numerical error: the numerics critic re-derived the Sorkin–Johnston reduction, the
measure, the UV tail and both IR prefactors analytically and reproduced the massless IR result with an independent
discretisation. Both rated the draft *sound with fixes*, and every fix was applied. The substantive ones: renormalisation
constants are a second free slot (§2, §5); the pre-registered N1 originally asked a question already answered by the
infrared divergence for ν ≥ 3/2 (§9); the horizon-state theorem's ν was conflated with the slow-roll mode index (§5);
Bombelli–Henson–Sorkin was applied to unembedded graphs (§3); "would describe a universe that is not ours" overstated a
fine-tuning debt (§3); "GRUT derives the initial state" asserted a derivation nobody has done (§7); the interference
argument was stated as established and its Gaussian law beyond weak coupling (§3, `v16`); the zero-temperature results of
the codex line were ignored by the original diagnosis (§0, §2).

## 9. What is not claimed, and the next pre-registered computations

Not claimed: a theory of everything, a new law, or a prediction. This document corrects the map and narrows the target.

1. **N1 — the horizon state for the massless field in quasi-de Sitter.** On the full field algebra a past-horizon
   Bunch–Davies-like state is IR-divergent for ν ≥ 3/2 (Ford–Parker 1977), so that question already has its answer. The open
   question is whether P1 selects, *without a free infrared scale*, either (a) a state on the shift-symmetric algebra
   generated by ∂φ, or (b) a unique infrared completion. *Success:* (a) or (b) exists, is Hadamard, and gives Bunch–Davies for
   the modes we observe. *Failure:* every such completion needs an infrared scale, and P1 cannot fix the state of the fields
   that seed structure. Either outcome is decisive for GRUT and is pre-registered here. Tools: the `v12` lemma and the
   expansion of `v14` around the Γ(−δ) pole.
2. **N2 — the second door.** Does any horizon principle in P0–P5 fix a renormalisation constant — Λ in particular —
   without a fit? None has been proposed on this branch; `06` (a₀ ∼ cH₀) is Milgrom's coincidence, not such a principle. A
   card in this direction must state the constant, the principle and the number before comparing with data.
3. **N3 — if the microstructure line continues,** it should switch to causal-order kinematics with Lorentzian amplitudes,
   pre-register the Lorentz-invariance benchmark (§3) and the entropy–action budget (`v16`), and compare with
   Mathur–Singh–Surya and Carlip–Carlip–Surya before computing anything new.
4. **N4 — Stage-3 Card 1** only for a rule that passes constraints 1–10 of `12`. No such rule exists today; the closest
   candidates are the subjects of N1 and N2.
5. Owed bookkeeping: the first-principles derivation of the zero-frequency invariant in `v13` 13.5; a literature check for
   r(ω) and the −(1 + 2 ln 2) coefficient.

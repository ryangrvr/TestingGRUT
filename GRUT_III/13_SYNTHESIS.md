# 13 — Synthesis of all the guidance: what is right, what cannot be right, and the one open door

**Status (2026-10-12).** Written at the owner's request to synthesise the guidance received so far, on the assumption that
some of it is wrong. The guidance: the GRUT III baseline and its Card-1 constraint sheet (`12`); the reviewer's fork
recommendation; ChatGPT's self-assessment of the program; and the codex graph line (fixed-bond entropy theorem, SCOUT-RG01,
SCOUT-RG02, CP-12 and its three audits, read from their branches). Method: every claim the synthesis rests on was checked
against primary sources by a six-group literature workflow and a codex-mathematics audit (§8 records verdicts and read
levels). Two new computations were added: `v15` (the Sorkin–Johnston state in a de Sitter slab) and `v16` (the
entropy–action budget). Where the checks contradicted my own earlier guidance, it is corrected here and in `09` and `12`.

## 0. The answer in one paragraph

Nothing in either line has found new physics, and the reason is now understood rather than felt. Both lines asked
*temperature* to do work it cannot do. With geometry given, thermality is a property of the state, so a thermal reading adds
nothing. Without geometry there is no bath, and Gibbs weights over a super-exponentially large configuration space cannot
select geometry; `v16` shows a local action's phase cannot either. The graph line's own target, bounded-degree spatial
graphs with a Lieb–Robinson cone in an external time, would if reached describe a universe with a preferred frame. Loop
effects would then make species' light cones differ by 0.1–10%, against measured equality to ~10⁻²⁰. What survives is
one open door, and it is narrow. A distinctive prediction can only come from a principle that *selects the quantum state*
without a temperature, because that is the one place where local covariance provably leaves a choice. Every known rule of
that kind fails or needs a scale the principles do not fix. The best parameter-free one, the Sorkin–Johnston state, is
computed here (`v15`): it violates GRUT's own P2 and contradicts the CMB over the inflationary slab. All of these failures
point to the same place, the infrared of massless fields in an accelerating universe. That is where the state is
genuinely underdetermined, where GRUT's horizon principle naturally speaks, and where observation is weakest.

## 1. What each source of guidance said, and its verdict

- **GRUT III baseline (`12`)**: a verified thermodynamic reading of GR + QFT, no distinctive prediction. *Stands.* Two of its
  lines were overstated and are corrected below (de Sitter KMS needs Hadamard; P2's state content is the Hadamard condition).
- **The reviewer**: bank the baseline and write the sixth postulate as Card 1 under the frozen charter; the options are not
  exclusive. *Right.* This synthesis turns "a sixth postulate that breaks KMS in a controlled way" into a precise target
  (§5) and adds three admissibility lines (`12` constraints 8–10).
- **ChatGPT's self-assessment**: much mathematics, no new law, progress over-reported; demand originality audits, constrained
  construction and an external benchmark. *Right in substance.* Its proposed replacement principle ("derive the configuration
  space from constraints") points the right way but is too vague to act on; §3 gives the precise requirement. The external
  benchmark it asked for exists and is decisive: Lorentz-invariance data (§3).
- **The codex graph line**: correct model-specific theorems; no spacetime. *Its mathematics stands (§8); its premise cannot be
  right* (§3).
- **My own earlier claim that "no state-selection rule has been proposed"** (`12` §3): *wrong.* Several exist (Sorkin–Johnston,
  states of low energy, horizon states, adiabatic vacua). Each buys uniqueness with a cost; §5 prices them.

## 2. The diagnosis: one wall, hit from both sides

The two lines look unrelated — a thermodynamic reading of continuum physics, and quantum graphs that should grow space.
They failed for the same reason.

- **GRUT III, with geometry given.** P1 assigns the bath the KMS temperature of the state's own modular flow wherever that
  flow is geometric. Quantum field theory already fixes every correlator of that state. So the thermal reading adds nothing:
  KMS by construction, no distinctive prediction (`12`; `v11`, `v14`).
- **The graph line, with geometry sought.** A Gibbs state at fixed temperature over a configuration space whose
  non-geometric part outnumbers the geometric part super-exponentially cannot select geometry with extensive energies, and —
  `v16` — a local action's phase cannot either. Thermal weighting cannot produce what the counting forbids.

In both, *temperature is being asked to do work it cannot do.* Once a geometry and a state are given, thermality is a
property of the state, not an extra input. Before there is a geometry, there is no bath — the universe has no outside —
so "a fixed positive temperature over pregeometries" is not a physical premise in the first place, and positive weights are
not what quantum gravity sums. Two corollaries follow, and they organise everything below:
1. Geometry has to come from **kinematics plus amplitudes**, not from thermal selection: a configuration space small enough
   (or an action large enough) that a Lorentzian amplitude can single out manifold-like structure.
2. A distinctive prediction has to come from **the state**, chosen by a principle that is not a temperature. That is the
   one slot in physics where local covariance provably leaves a choice open, and where GRUT's horizon relationality (P1)
   could be the deciding principle.

## 3. The graph line: what is right, what cannot be right

**Right, and worth keeping.** The codex theorems are correct as stated (independently re-derived in the workflow, §8):
conserving the number of connections removes the dense-edge failure, but at fixed positive temperature a family of lattice
templates (at most N!·e^{O(N)} labelled members) has probability at most exp[−(ν−1)N ln N + O(N)] among graphs with
Q = νN connections, even with the full non-commuting gauge–fermion–rewiring Hamiltonian (relative-entropy bound). The CP-12
back-reaction (a charged connection is removed at order T⁴ instead of T²), the Gauss-law distance diagnostic, the
O(T²/N) operational window and the coherent no-assembly bound are sound model-specific results.

**Rediscovered, not discovered.** The central obstruction is thirty years old in other notation. Chen & Plotkin (Phys. Rev.
D 87, 084011, 2013) found numerically that the transition temperature of fixed-edge graph models goes to zero as N grows;
Kleitman & Rothschild (1975) is the same obstruction for causal orders; Euclidean dynamical triangulations met it as the
crumpled and branched-polymer phases. Quantum graphity (Konopka, Markopoulou & Severini 2008; Hamma et al. 2010) met the
locality and equilibrium problems the RG02 and CP-12 audits found. Trugenberger and Kelly–Trugenberger–Biancalana obtained
geometric phases only by restricting the configuration space by hand (regular bipartite graphs) and rescaling the coupling
with N — the two moves the codex charter forbids, which is why the codex line cannot reach the same conclusion.

**The precise reason (v16).** Selection needs an *entropy–action budget*. A weight can make a geometric class dominate only
if what it can vary across the competing class is at least that class's excess log-count.
- Boltzmann weights with extensive energies (|E| ≤ CN): the codex theorem — excess (α−1)N ln N wins.
- Lorentzian weights e^{iS} with a local action on sparse graphs: S has O(N) roughly independent terms, so the phase
  suppresses the class sum by about exp(−λ²Var S/2) = exp(−O(N)) (v16 checks the Gaussian law to 10% and Var S ∝ N).
  Still loses to N ln N. Interference is not a loophole for local actions.
- The two ways out that are known: a configuration space whose counts grow only exponentially (triangulations of fixed
  topology: proven in 2D, for locally constructible 3D triangulations, and for 3D causal triangulations; open in general
  3D and in 4D), or a nonlocal action with ~N² terms against N² entropy (causal-set theory: the dominant layered orders are
  suppressed by the link or BDG action — Loomis–Carlip 2018 for two-level orders; Mathur–Singh–Surya 2021 and
  Carlip–Carlip–Surya 2023–24 for Kleitman–Rothschild/layered orders — at leading order and in periodic parameter windows;
  no paper yet compares these sums with the manifold-like contribution, and the full two-level class is not exponentially
  small, so this is necessary-and-promising, not solved).
- One sharpening of the codex result (graph-model check, §8): if vertex relabellings are treated as gauge (unlabelled
  counting), the competing class is super-exponential for *every* ν > 0 — at least exp(νN ln N − 2νN ln ln N + O(N)) —
  so CP-12's "ν = 1 is open" closes in that reading. In labelled counting it stays open.

**What cannot be right: the target itself.** CP-12's stated goal is a bounded-degree, spatially organised graph on which
matter propagates inside a Lieb–Robinson cone in external time. Three independent facts say that target, if reached, would
describe a universe that is not ours:
1. *Bounded degree picks a frame.* Bombelli, Henson & Sorkin (Mod. Phys. Lett. A 24, 2579, 2009) prove there is no
   Lorentz-equivariant measurable rule assigning a preferred direction to a Poisson sprinkling of Minkowski space; as a
   corollary no finite-valency neighbour graph can be attached equivariantly. A Lorentz-invariant discreteness has to be
   nonlocal: in a sprinkling every element has infinitely many links. A Lieb–Robinson cone in an external time is a cone
   in the lattice frame by construction.
2. *A frame at the Planck scale does not stay there.* Collins, Perez, Sudarsky, Urrutia & Vucetich (PRL 93, 191301, 2004):
   a Lorentz-violating ultraviolet regulator feeds, through ordinary loops, unsuppressed dimension-4 terms that shift each
   species' limiting speed by c → c(1 + ξ/4), ξ = (g²/6π²)[1 + …] — 0.1% to 10% for Standard Model couplings, independent of
   the cutoff scale. No custodial symmetry is known except low-scale supersymmetry (not observed); logarithmic infrared
   running (Chadha–Nielsen) is far too slow.
3. *Experiment.* Limiting-speed differences between species are bounded at ~10⁻²⁰–10⁻²¹ (electron–photon, proton–photon,
   one-sided astrophysical; Kostelecký–Russell tables, 2026 edition), photon birefringence at ~10⁻³⁵, and linear
   energy-dependent photon speeds above ~2 E_Planck (Vasileiou et al. 2013, conservative; Fermi GRB 090510 > 1.2 E_Pl).
   A bounded-degree pregeometry therefore owes an explanation of 10⁻²⁰ equality of species' light cones. That is the
   external, unfitted benchmark the codex line was asked to find, and it is decisive.

**External time.** Every codex model (CP-11, CP-12, RG01, RG02) evolves in an external parameter t with a Hamiltonian
bounded below (structural survey, §8). General relativity has no such time: its Hamiltonian is a constraint. A model whose
causal cone is a Lieb–Robinson bound in an external clock has supplied the preferred foliation that general covariance
forbids. This is the same obstruction as (1), seen from the dynamics.

**Verdict on the line.** Its obstructions are real and, now, understood. Its premise — undirected adjacency, Gibbs
weights, external time, bounded degree — has each of its four ingredients independently excluded or at best unexplained.
Continuing to add mechanisms inside that premise will keep producing correct theorems about models that cannot be our
universe. If the line continues, it should change kinematics: a causal order (directed, acyclic) plus a counting measure
— which by Hawking–King–McCarthy/Malament fixes a Lorentzian metric up to its conformal factor in d > 2 — with a
Lorentzian amplitude and an action whose budget matches the entropy. That is causal-set theory's choice, made for exactly
these reasons; GRUT should build on it rather than re-derive its obstacles.

## 4. Corrections to the GRUT III guidance (mine)

- **De Sitter KMS needs Hadamard.** `12` said the de Sitter vacuum is KMS "by de Sitter invariance and positivity". De Sitter
  invariance alone allows the Mottola–Allen α-vacua, which are not thermal for geodesic observers. Bunch–Davies is the unique
  de Sitter-invariant pure Gaussian *Hadamard* state for m² > 0, and for the massless minimally coupled field no de
  Sitter-invariant state exists at all (Allen 1985; Hollands–Wald 2015). Constraint 1 of `12` is corrected.
- **P2's state content is the Hadamard condition, as a scaling-limit statement.** Jacobson's 1995 "local equilibrium" has a
  geometric part (an instantaneously stationary local Rindler horizon) and a state part ("the usual short-distance vacuum
  fluctuations"). Only the state part is the Hadamard condition. Exact KMS at κ/2π follows from Hadamard only at bifurcate
  Killing horizons (Kay–Wald 1991); at other horizons it holds in a near-horizon scaling limit (Kurpicz–Pinamonti–Verch 2021;
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

## 5. The one slot a principle could fill — and the test of the best existing candidate (v15)

Every distinctive claim has to change the state, the operator content or the identification (`12` §3). Changing the
operator content is conventional physics with a price. Changing the identification abandons GRUT III. That leaves the
state, and there the field has a precise, known gap: **no rule selects a physical state of a quantum field in a general
curved spacetime that is at once locally covariant, Hadamard and free of arbitrary choices** (Fewster–Verch 2012, Thm 6.13; §4). In de Sitter the Bunch–Davies state is selected by symmetry plus Hadamard; in FRW every known rule carries a
choice — an initial time (adiabatic vacua), a smearing function (Olbermann's states of low energy), a smoothing
function (Brum–Fredenhagen 2014) — the same kind of hidden scale that `v13` showed can move a first-order coefficient
anywhere. A GRUT principle that selected the state with none of these would be a genuine addition, and this is the place
P1 (the bath is the observer's causal complement) naturally speaks, because a selection rule that local covariance cannot
supply has to use global data.

**The best existing candidate, computed.** The Sorkin–Johnston (SJ) prescription (Johnston 2009; Afshordi, Aslanbeigi & Sorkin 2012) is parameter-free and
temperature-free: keep the positive part of the Pauli–Jordan operator on the region's L² space. `v15` reduces it exactly,
mode by mode, in a de Sitter slab ending at a_e H (two integrals per mode; the SJ condition is ∫v_k² a⁴dη = 0), checks the
reduction against the positive part of the discretised operator to 4e-12, and finds:
- **Ultraviolet:** |β_k|² → (a_eH)²/(16k²) for both the conformal and the massless field — a power law, so the state is not
  Hadamard (Fewster–Verch 2012 proved the same for ultrastatic slabs), and the |β|² part of its energy density grows like Λ² with the cutoff.
  It violates P2.
- **Infrared:** modes outside the horizon at the region's end are suppressed, R_k = |v_k|²/|u_k|² = √(π/2)(k/a_eH)^{1/2}
  (conformal) and √(π/2)(k/a_eH)^{3/2} (massless). Over the inflationary slab every observable mode would have its power
  multiplied by (k/a_eH)^{3/2}: a tilt shift of +1.5 against Planck's n_s − 1 = −0.035 ± 0.004. Excluded.
- **Teleology:** extending the region's future boundary from a_eH = 1 to 100 changes the state of a k = 10⁻³ mode at a fixed
  earlier time by three orders of magnitude. The state today would depend on where "the region" is declared to end.

So the only known rule of the right kind fails GRUT III's own P2 in a sharp region and contradicts the CMB over the
inflationary slab. Smoothing it restores Hadamard and reintroduces a free function. This is the honest state of the art,
and it sharpens what Card 1 has to be: a rule that is Hadamard, causal (no dependence on a future boundary), free of
arbitrary functions, and that reproduces Bunch–Davies for the modes we observe.

**Where this sits in the literature.** Sorkin–Johnston states are not Hadamard on ultrastatic slabs for almost every slab
width, and are not locally covariant (Fewster–Verch 2012). In de Sitter the limiting Sorkin–Johnston state is always an
α-vacuum (Aslanbeigi–Buck 2013): Bunch–Davies on the global patch in 4D only for m ≥ 3/(2ℓ); for lighter fields a
non-Hadamard α-vacuum, or ill-defined (including the conformally coupled massless field); undefined on the Poincaré patch
below that mass. The finite-slab results of `v15` for the light fields are the finite-region face of the same facts; the
explicit UV tail and IR exponents are direct computations, not quoted. Afshordi, Aslanbeigi & Sorkin (2012) found that in a
radiation-era FRW the massless Sorkin–Johnston state is the conformal vacuum, and they stress the region dependence.

**The other known rules, priced.** States of low energy (Olbermann 2007; Them–Brum 2013) are Hadamard for every choice but
depend on a smearing function f. Smoothed Sorkin–Johnston states (Brum–Fredenhagen 2014) are Hadamard but depend on a
smoothing function. Adiabatic vacua depend on an initial time and an order. The horizon state of Dappiaggi, Moretti &
Pinamonti (2009) comes closest to GRUT's P1. It is built from data on the past cosmological horizon, is unique and invariant
under the horizon's symmetry group, equals Bunch–Davies in exact de Sitter, and is proven Hadamard for spatially flat FRW
that is asymptotically de Sitter in the past. **But the proof needs Re ν < 3/2** (or Re ν < 1/2 under weaker asymptotics)
and excludes the massless minimally coupled field, ν = 3/2. In quasi-de Sitter that field has ν = 3/2 + ε/(1−ε) > 3/2 (the
`v12` lemma), outside the theorem. So the one GRUT-native rule that works is not yet established for the fields that seed
cosmic structure.

## 6. The unifying finding: every route ends at an infrared scale

Every attempt in this program to get beyond the conventional answer has stopped at the same kind of object: a scale the
principles do not fix, in the infrared of massless fields in an accelerating spacetime.

| Route | Where it stopped |
|---|---|
| `v13` zero-mode "temperature" | the coarse-graining scale σ moves the O(ε) coefficient anywhere |
| `v14` massless sector | the IR constant and the secular term, which the Wigner convention puts at ω = 0 |
| `10` secular-IR survey | the graviton dispute is over which IR quantities are observable |
| Sorkin–Johnston (`v15`) | the region's future boundary, through k/a_eH, decides the superhorizon state |
| low-energy and smoothed states | a smearing or smoothing function |
| horizon state (DMP) | the Hadamard theorem stops at ν = 3/2; the massless field in quasi-de Sitter is beyond it |

This is not a coincidence. For ν ≥ 3/2 the Bunch–Davies two-point function has the Γ(−δ) pole that `v14` expanded around,
the infrared divergence of a massless field in an accelerating universe. The state of those modes is fixed by when and how
the acceleration began, not by local physics. So the infrared is where conventional theory genuinely leaves the state
open. It is where a horizon-relational principle like P1 could act. And it is the one place observation cannot reach well:
modes that left the horizon near the start of inflation are today at or beyond the largest observable scales, where cosmic
variance dominates. Planck finds every non-Bunch–Davies bispectrum template consistent with zero, and any oscillatory
modulation of the power spectrum below about 3%.

**Consequence.** A GRUT prediction, if one exists, lives in the lowest multipoles of the CMB and the largest scales of
structure. It must be made with its scale fixed by the principle before looking, not fitted to the ~2σ low-ℓ deficit. With
the principles as they stand, the scale is not fixed, so there is no prediction.

## 7. What a theory built on these lessons would look like

Three layers, each with the benchmark it must meet. None is complete; this is the map, not the territory.

1. **Microstructure: causal order plus counting, not adjacency.** In d > 2 a causal order fixes a Lorentzian metric up to a
   conformal factor (Hawking–King–McCarthy 1976; Malament 1977), and counting supplies the volume ("order + number =
   geometry"). A Poisson sprinkling is statistically Lorentz invariant, so Lorentz invariance is built in rather than owed.
   The price is nonlocality: infinitely many links per element. Weights must be Lorentzian amplitudes with an action whose
   budget matches the entropy. *Benchmark:* manifold-like orders must dominate the layered classes. Today the dominant
   layered classes are suppressed at leading order in periodic parameter windows, and the comparison with manifold-like
   orders is open.
2. **Low-energy dynamics: general relativity from P2.** Jacobson 1995 and the entanglement-equilibrium version (Jacobson 2016)
   give the Einstein equation from local horizon thermodynamics, given a Hadamard-type state and a finite universal entropy
   density. *Benchmark:* GR's tests, which GRUT III already passes because it is GR.
3. **The state: Hadamard (P2), selected horizon-relationally (P1).** Where a past cosmological horizon exists and the theorem
   applies, this is the DMP state: Bunch–Davies in de Sitter. GRUT's principle then *derives* the standard initial state of
   inflationary cosmology instead of assuming it. That is explanatory, not predictive. *Benchmark:* Planck's
   non-Bunch–Davies and oscillation bounds, and n_s.

What GRUT adds, honestly stated: P1 as the selection rule in layer 3, and an interpretation of layer 2. It adds nothing yet
to layer 1, where causal-set theory has the right kinematics and GRUT has no dynamics of its own.

## 8. Verification record

CODEX_PLACEHOLDER

## 9. What is not claimed, and the next pre-registered computations

Not claimed: a theory of everything, a new law, or a prediction. This document corrects the map and narrows the target.

1. **N1 — the horizon state for the massless field in quasi-de Sitter.** Does a past-horizon (DMP-type) state exist for
   ν = 3/2 + ε/(1−ε), and is it Hadamard, without an infrared scale? *Success:* it exists, is Hadamard, and gives
   Bunch–Davies for the modes we observe; GRUT's P1 then derives the inflationary initial state. *Failure:* it needs an
   infrared scale, and P1 cannot fix the state of the fields that seed structure. Either outcome is decisive for GRUT and is
   pre-registered here. Tools: the `v12` lemma and the expansion of `v14` around the Γ(−δ) pole.
2. **N2 — if the microstructure line continues,** it should switch to causal-order kinematics with Lorentzian amplitudes,
   pre-register the Lorentz-invariance benchmark (§3) and the entropy–action budget (`v16`), and compare with
   Mathur–Singh–Surya and Carlip–Carlip–Surya before computing anything new.
3. **N3 — Stage-3 Card 1** only for a rule that passes constraints 1–10 of `12`. No such rule exists today; the closest
   candidate is the subject of N1.
4. Owed bookkeeping: the first-principles derivation of the zero-frequency invariant in `v13` 13.5; a literature check for
   r(ω) and the −(1 + 2 ln 2) coefficient.

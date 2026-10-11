# 01 — Principles of GRUT III (revision 2)

Tiers: **POSTULATE** (assumed; counted), **IMPORTED** (an established theorem or result brought in; zero credit to
GRUT), **DERIVED** (computed in `calc/` from the postulates plus imports), **DATUM** (measured constant; not
counted), **CONJECTURE**. Revision 2 corrects revision 1 in two ways flagged by the reviews: imports are labelled
as imports, and P1 is stated in the algebraic form that type III local algebras require.

## Data (taken from nature)

ħ, c, G, k_B; the cosmological constant Λ; the particle spectrum (masses, spins, multiplicities, gauge group).
G and ħ are data by the charter's own ruling on Λ; **this is also where the gravitational coupling enters the
framework** (through P2 and P3 below). GRUT III does not derive G and says so (05 O1).

## P0 — Covariant quantum field theory (POSTULATE, +1)

Physics is a local relativistic quantum field theory on a Lorentzian manifold: a net of von Neumann algebras A(O)
obeying Einstein causality and covariance, with a vacuum satisfying the spectrum condition and states of Hadamard
form. **This single postulate bundles a great deal**: Lorentzian signature, causal locality, the Hilbert-space
structure, and the field content (the data). It is one assumption in the sense that it is one framework; it is not
"less" than the record's rung-1 package, it is the same package named honestly. It replaces rung 1 (+4) and rung 6's
"quantization" (+1) in the old ledger without earning anything.

## P1 — Horizon relationality (IDENTIFICATION, 0; stated algebraically)

For an observer whose accessible region is the causally complete O, the **medium** is the triple
(A(O), ω|_{A(O)}, σ^ω_t): the local algebra, the restriction of the global state to it, and the Tomita–Takesaki
modular flow of that restricted state. The **bath** is the commutant A(O)′ (= A(O′) where Haag duality holds).

What P1 does *not* do (R8): local algebras are type III₁, so there is no tensor factorisation H = H_O ⊗ H_{O′},
no reduced density matrix, and P1 supplies **no influence functional and no dynamics**. What it does: by
Tomita–Takesaki the restricted state is KMS at β = 1 in the modular parameter automatically, so "the medium is
thermal" is a tautology; the *physical* content is whether the modular flow is **geometric**, which is exactly when
a temperature is a number:

- Minkowski wedge: modular flow = boosts, β = 2π/a (Bisognano–Wichmann 1975/76) — **IMPORTED**.
- Bifurcate Killing horizon: unique invariant Hadamard state, KMS at β = 2π/κ (Kay–Wald 1991) — **IMPORTED**.
- De Sitter, geodesic observer: Bunch–Davies is thermal at T = H/2π (Gibbons–Hawking 1977) — **IMPORTED**;
  worldline correlators in closed form in `v1` (D1–D2).
- Double cones in a CFT (Hislop–Longo) — **IMPORTED**; generic regions: non-geometric, no temperature.

The influence functional for a probe coupled through P3 is then *computed* (Feynman–Vernon), and the record's
Gaussian in-in action S_IF is its quadratic approximation. P1 replaces the record's `rung1_ontology_finite_memory`
and `background_time_translation_flow` by naming the same objects; it earns nothing by itself.

## P2 — Local equilibrium of the medium (POSTULATE, +1; G enters here)

Across every local causal horizon the medium's entropy is S = ηA and heat flux obeys Clausius δQ = T δS with T
the KMS temperature of P1. Area-proportionality of the entanglement entropy is **IMPORTED** (Bombelli–Koul–Lee–
Sorkin 1986; Srednicki 1993) with a UV-divergent coefficient; **η = 1/(4Għ) is inserted** (the datum G) — its
finite value is the open problem O1, and `v9` shows that one-loop induced gravity with Standard-Model content does
*not* supply it. The Clausius relation is the one structural postulate. Consequence (IMPORTED derivation, Jacobson
1995, reproduced in `v2`): G_ab + Λg_ab = 8πG T_ab with Λ an integration constant. GR is **recovered**, not derived.

## P3 — Universal coupling (POSTULATE, +1; κ enters here)

The medium couples to matter only through T_μν: the vertex is (κ/2) h_μν T^μν with κ² = 32πG. The response
kernels are then the ⟨TT⟩ functions and nothing else; the absorptive kernel Im Π(s) = c_eff(s) G s² is **DERIVED**
from this vertex plus the spectrum (`v4`), as are the tensor/scalar parts of the kernel (no `p_tt` ansatz). The
coupling itself is inserted; its universality is the operational equivalence principle.

## P4 — Time orientation (POSTULATE, +1)

The Hadamard condition is imposed with positive frequencies toward the future: the traced-out horizon is the
future horizon, β > 0. This is the record's Past Hypothesis, located; it is not derived.

## P5 — Outcome statistics (POSTULATE, +1)

Outcome probabilities depend only on the effect measured (non-contextuality); Gleason (dim ≥ 3) and Busch (POVMs)
then give Tr(ρE) — **IMPORTED**. Revision 1's suggestion that P5 might follow from P0 + P1 via
Hughston–Jozsa–Wootters is **withdrawn** (R5): HJW presupposes the measurement framework. The reconstruction
literature (Gleason; Zurek's envariance; Deutsch–Wallace; Masanes–Galley–Müller 2019, who derive POVMs, Born and
state update from unitary quantum mechanics plus operational closure, associativity and finite state estimation —
**without** assuming non-contextuality, i.e. a stronger result than P5's route) trades the postulate for other
assumptions; GRUT III adds nothing here.

## Accounting

| | |
|---|---|
| Framework postulates | P0, P2, P3, P4, P5 = **5** |
| Identifications (no content by themselves) | P1 |
| Imported theorems (zero credit) | Bisognano–Wichmann; Kay–Wald; Gibbons–Hawking; Hislop–Longo; entanglement area law; Jacobson EOS; Eling–Guedens–Jacobson membrane transport; Huygens for conformal fields; Gleason/Busch; Allen–Folacci dS propagator; Starobinsky–Yokoyama stochastic diffusion |
| Derived in `calc/` | the closed-form dS worldline spectra; the species kernel; the smeared-response theorems; the gravitational-sector random walk; the quasi-dS KMS defect; the induced-gravity matching |
| Data | ħ, c, G, k_B, Λ, spectrum |

The record's +16 is not "compressed" to +5. Seven of its nodes are now imported theorems (zero credit), two are
data, two are identifications; five framework postulates remain. **No kernel is supplied anywhere; nothing is
earned by that fact alone**, because the kernels that replace it are standard physics.

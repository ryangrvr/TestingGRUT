# CP10: derive the reference's actual geometry feedback

2026-10-11. Additive analysis of frozen CP10, parent abf99c91b7560602caf089ed483079f0fde7dc2f. No new generating law, no CP11, no Stage3 card. Fixed-geometry controls below are explicitly not solutions of the full fluctuating-geometry Hamiltonian. They isolate its supplied reference mechanism. No fundamental-discovery or literature-priority claim.

## Exact full-model facts

Write H10 = H09 ⊗ I + HR + gW, with W = Σ_i (I-Z_i)⊗|i><i| and HR the configuration-controlled graph Laplacian. Let |u> = N^(-1/2)Σ_i|i>.

Every ordinary graph Laplacian annihilates |u>, including disconnected graphs. Consequently HR(|ψ>⊗|u>)=0 for EVERY geometry–matter state |ψ>. At g=0, the entire CP09 spectrum has an exact copy in CP10; E0,10=E0,09. The added propagation term alone supplies no ground-state selection. Raising k cannot remove this common zero mode.

Since 0≤W≤2I, variational positivity gives the uniform bound

    E0,09 ≤ E0,10 ≤ E0,09 + 2g.

At frozen g=1/4, this energy correction is at most 1/2 at every N. This does NOT bound the ground-state probability distribution: arbitrarily small perturbations can rearrange nearly degenerate ground states. In particular it is not a thermodynamic no-selection theorem. It explains why one finite reference is not automatically an extensive back-reaction mechanism.

## A solvable conditional energy surface

Hold geometry fixed and isolate the reference and onsite matter Hamiltonian:

    h_G = Ω Σ_i(I-X_i) + k I_m⊗L_G + g Σ_i(I-Z_i)⊗|i><i|.

This omits CP09 configuration hopping and its degree diagonal. Do not add its ground energies to CP09 as an exact Born–Oppenheimer potential without a separation-of-scales argument.

If G' adds an edge (i,j),

    h_G' - h_G = k I_m⊗(|i>-|j>)(<i|-<j|) ≥ 0.

Thus its lowest energy is monotonically NONDECREASING under edge addition for any g≥0. This is nonperturbative. The reference sector energetically prefers the empty graph over every other graph; it cannot alone favor a connected physical geometry over emptiness. If transitivity demands several added edges the same conclusion holds term by term.

For the empty graph, the occupied-site block has one spin with Hamiltonian Ω(I-X)+g(I-Z); other spins are in X+. Hence exactly

    E_empty = Ω+g-sqrt(Ω²+g²).

A uniform reference and all-X+ trial state has energy g on every graph. Therefore

    Ω+g-sqrt(Ω²+g²) ≤ E_G ≤ g.

This gives an explicit, N-independent range for the entire conditional geometry energy surface.

## Derived leading selection coefficient

For connected G, at g=0 the reference ground state |u> is unique and the matter ground state is |+>^N. The first-order correction is g. Acting with the remaining perturbation creates one X-flipped spin and a localized reference amplitude:

    (W-I)|+>^N|u> = -1/sqrt(N) Σ_i Z_i|+>^N |i>.

Different spin flips are orthogonal and cost 2Ω. Resolving reference excitations gives

    E_G(g) = g - (g²/N) Tr[(kL_G+2Ω I)^(-1)] + O(g⁴).

The absence of the cubic term follows because conjugation by Π_i X_i maps h_G(g)-gI to h_G(-g)+gI. Its isolated ground eigenvalue is even in g after subtracting g. The expansion is at fixed finite connected graph; no uniform large-N radius of convergence is claimed. Disconnected graphs have a degenerate reference vacuum and require degenerate perturbation theory.

This is an explicit microscopic spectral response: geometry enters through the full graph Laplacian resolvent. It rewards soft reference modes, but softness does not distinguish spacetime from stars, bottlenecks or disconnected configurations.

At Ω=k=g=1/4, exact conditional ground energies are:

| N | Empty | Height-two star | Balanced height two | Total causal order |
|---|---:|---:|---:|---:|
| 2 | .146446609 | .166938375 | .166938375 | .166938375 |
| 3 | .146446609 | .172471073 | .172471073 | .181037978 |
| 4 | .146446609 | .174552015 | .182638838 | .191223228 |
| 5 | .146446609 | .175464834 | .187312211 | .198878016 |

All listed graphs are comparability graphs of actual posets. The star is a valid height-two order. The total causal order has COMPLETE comparability graph, not a nearest-neighbor chain.

## A propagation failure visible without diagonalizing the full model

For a total order, L = N I - 11^T. The reference-only propagator is

    U(t) = |u><u| + exp(-ikNt)(I-|u><u|).

For every distinct pair i,j,

    |<j|U(t)|i>|² = 4 sin²(kNt/2)/N².

There is no dependence on their separation in the total causal order. Ordering more intermediate events between source and detector does not delay propagation. This is a test of the supplied hopping term, with geometry held fixed and matter coupling absent; it is not the complete interacting detector amplitude. Even its total departure probability is ≤4(N-1)/N², illustrating why rapid frequencies alone do not imply useful transport.

Graph-Laplacian quantum walks and complete-graph transport are established mathematics; see Razzoli, Paris and Bordone, *Universality of the fully connected vertex in Laplacian continuous-time quantum walk problems*, https://arxiv.org/abs/2202.13824. The CP10 application and perturbation calculation here are author-level deductions, not claims of new universal physics.

## Execution and next construction requirements

The committed check.py independently assembles h_G, checks exact empty energies and variational bounds, compares exact diagonalization with the derived second derivative at g=.004,.002,.001, and compares matrix exponentiation with the complete-graph transition formula, for N=2..5. All assertions passed; outputs are results.json. Small-g points estimate a derivative of the frozen operator family; they are not alternate candidate parameter selections. Frozen g=.25 energies are separately listed.

The scientific consequence is a concrete design constraint: merely adding more copies of this reference can magnify a preference for low-connectivity graphs without selecting spacetime. Replacing all comparabilities by covers would introduce a different microscopic law and still would not derive Lorentz locality. Neither is an earned repair.

A next candidate must supply a physical competition that stabilizes a connected phase while preserving an operational propagation limit. It must declare what supplies connectivity, how the propagation kernel handles causal-order nonlocality, and whether the same interactions generate the geometry energy. A graph chosen because its modes look spatial would assume the sought result. The present calculation identifies the missing competition quantitatively rather than claiming that a responsive reference has already supplied it.

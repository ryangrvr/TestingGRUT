# CP-12: local gauge capacity and the cost of preparing it

**2026-10-14. Frozen law unchanged. Additive continuation of the holonomy audit.**
Parent: `ae91c9beb0d8dda0f7692ea8bab0c5f808355530`.

Higher connection density does provide extensive local circulation degrees of
freedom on prepared spatial hosts. I tested this explicitly at Q=2N and Q=3N,
and checked their leading Wilson-loop response with the full frozen Hamiltonian.
But increasing Q does not by itself make those hosts typical. A new conditional
bound covers a much broader geometric family than the earlier template-counting
argument, and it also constrains a class of nonequilibrium preparations.

## 1. The supplied audit reproduces

The exact committed `cp12_holonomy_locality_gate/check.py` was executed verbatim.
It passes the unicyclic graph counts and Gauss-incidence controls. The theorem
that a connected Q=N graph has one independent circulation is correct. The
parent's small rank calculation uses floating elimination; the new higher-density
rank controls below use exact rational elimination instead.

The stronger interpretation is **local circulation capacity**, not connectivity
alone. This agrees with the audit's qualification that the Q=N ground state
cannot directly be read as a conventional three-dimensional lattice gauge system.
The reported small-system connectedness and longer-cycle bias remain intact.
I have not independently re-proved the earlier near-degree-two ground-state
exclusion; the audit's bridge corollary retains that stated dependency.

## 2. A positive test at higher density

| Prepared host | N | Q | All independent cycles | Independent elementary-square circulations |
|---|---:|---:|---:|---:|
| Periodic 6 x 6 grid | 36 | 72 | 37 | 35 |
| Periodic 5 x 5 x 5 grid | 125 | 375 | 251 | 248 |

The residual 2 and 3 circulations wind around the tori. These counts are exact
rational ranks of explicit integer loop vectors. Coordinates construct the
benchmark hosts; they do not enter or change CP-12's Hamiltonian.

For the normalized Hermitian Wilson probe X_c=(W_c+W_c^dagger)/sqrt(2), a neutral
sea on a supplied cycle c has the full-law retarded response

\[
\boxed{i\langle[X_c(T),X_c]\rangle
=2|c|\left[\epsilon_E-\frac{\Gamma(M-Q)}M\right]T+O_N(T^3).}
\]

This includes the reduction of the rewiring diagonal when flux locks the cycle.
The two hosts yield exact slopes **226/175** and **192/155**. The computation
streams every one-action output, including 2,765,625 rewiring destinations for
the 125-vertex neutral state. It is not a many-body eigensolve. A separate complete
1,875-state matrix and finite source impulses recover the same commutator identity.

This is a genuine internal U(1)-gauge response. It supplies neither a photon
pole nor a propagation velocity; a chosen cycle and host are still supplied.

## 3. New selection theorem: many local loops remain thermally rare

Fix R, and let K_R(G) count the largest collection of vertex-disjoint simple
cycles of length at most R. This is a deliberately broad target: **all** graphs
with K_R>=alpha N are included, whether spatial, distorted, fragmented or otherwise.
No lattice-template count is assumed.

An elementary fixed-Q counting bound is

\[
U_{N,Q}(K_R\ge k)\le\frac{A_{N,Q,R}^k}{k!},\qquad
A_{N,Q,R}=\sum_{r=3}^{R}\frac{(NQ/M)^r}{2r}.
\]

At Q=nu N with fixed nu and k=ceil(alpha N), its logarithm is at most
-alpha N log N+O(N). Exact enumeration of 127,972 graphs checks the packing
expectation and union-bound steps. Even-coordinate square and cubic tori have
N/4 disjoint elementary squares. More generally, bounded degree plus a positive
density of vertices on bounded-length cycles implies an extensive packing.

For the **complete quantum Gibbs state**, including unbounded integer electric
fields and coherent rewiring, I obtain

\[
\boxed{P_\beta(K_R\ge\alpha N)=O(1/\log N)\longrightarrow0}
\]

at every fixed positive temperature in the full labelled Hilbert space.
This is proved using a diagonal electric reference Gibbs state and quantum
relative entropy; the proof never replaces the full Gibbs state by a classical
graph distribution. It also holds when beta_N grows more slowly than log N,
provided beta_N is bounded away from zero.

This closes part of the earlier “perhaps the spatial family is much broader”
loophole: the target is now defined by many local circulation paths, rather than
a chosen list of geometric templates. It still does not cover all possible
operational or coarse-grained notions of space.

## 4. A new preparation constraint, valid at every time

Let sigma_beta be the electric-reference Gibbs state used in the proof. For
closed evolution under the frozen full Hamiltonian,

\[
D(\rho_T\|\sigma_\beta)
=D(\rho_0\|\sigma_\beta)+\beta(\langle V\rangle_0-\langle V\rangle_T),
\]

where V is hopping plus rewiring and its spectral width is O(N). Therefore,
if D(rho_0||sigma_beta)=O(N),

\[
\boxed{\sup_T P_T(K_R\ge\alpha N)=O(1/\log N).}
\]

The time supremum includes times increasing with N. In particular, this excludes
assembly of that loop-rich class from the **incoherent equal mixture of all
neutral fixed-Q graphs**. Waiting longer does not evade the result for this
preparation. The theorem assumes closed evolution without postselection.

The corresponding **coherent pure** equal superposition escapes the bound:
its initial relative entropy exceeds that of the mixture by exactly

\[
\log\binom MQ=\nu N\log N+O(N).
\]

Both preparations have extensive full energy, although the coherent state's
rewiring energy is lower. They have identical initial graph probabilities but
different quantum information resources. The pure preparation is not thereby
shown to assemble geometry. Other low-entropy preparations can also escape;
coherence is one option, not the uniquely required resource.

For a nonvanishing final probability p_* of the loop-rich class, the proof gives
the necessary initial budget D(rho_0||sigma_beta)>=p_* alpha N log N-O(N).
This is a concrete requirement on a proposed preparation, not an invented term
in the Hamiltonian.

## 5. Checks and remaining scope

`verify.py` freshly executes four checkers: the unchanged parent audit, exact
new graph controls, exact higher-density full-law moments, and complete small
quantum Gibbs/quench/response calculations. It does not merely inspect saved
PASS labels. See `VERIFICATION.json`, the three result files and `MANIFEST.json`.

The small quantum matrices have dimensions 688 and 1,875 and noncommute with the
tested graph projectors. Gibbs identities, binary data processing, interaction
bounds, and the all-time relative-entropy conservation identity pass at three
reference temperatures and four evolution times. These finite-cutoff checks
validate algebraic steps. The infinite-rotor and large-N conclusions rest on the
proof in [DERIVATION.md](DERIVATION.md), not numerical extrapolation.

The new bounds do **not** decide the Q=2N or Q=3N ground-state phase. They do not
apply automatically to a permutation-gauge quotient, external entropy removal,
loops whose lengths grow with N, or a geometry encoded only in coarse-grained
quantum correlations. No photon, graviton, relativity or dimension-selection
claim is made. The constants are too loose for useful small-N probability bounds.

The next justified constructive test is a specified low-entropy preparation in
one density sector, followed by the full coupled evolution, measuring local loop
packing, intrinsic volume growth and relational propagation together. The new
result identifies a necessary resource for that test and demonstrates actual
gauge capacity at higher density. It does not supply the missing spatial phase.

The entropy and data-processing methods are established physics/mathematics;
primary references and all assumptions are listed in DERIVATION.md. Independent
partner audit of the new proof remains outstanding.

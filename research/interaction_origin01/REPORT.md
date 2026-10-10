# What the GRUT substrate does not fix about interactions

**10 October 2026 · bounded author-level standard-model audit · independent review pending.**

**Result:** fixing the complete vacuum/one-particle channel of a chosen quantum
lift does not fix its interacting completion. A finite, positive,
number-conserving, two-site family has identical lower-sector dynamics for every
time and every intervention confined to those sectors, while its measurable
two-particle addition frequency varies continuously. Explicit energy-conserving
finite collisions realize the same lower-sector equivalence. This is a
conventional countermodel to a derivation from restricted substrate data, **not
a new law, an unexplained freedom of complete dynamics, or a candidate selector**.

Recon01–04 are retained as completed author-level conventional investigations.
The owner has closed record allocation as the principal discovery route. This
packet does not extend those record experiments or amend their frozen files.

## 1. What the repository already settled

The source files are pinned to starting commit
`ab23f21cb6af9dd47483780fbfc6843680b44101`; hashes and read scopes are in
`SOURCES.json`.

| Source | Binding boundary used here |
|---|---|
| `S5_GENERATOR_ORIGIN_01.md`, `S5_GENERATOR_ORIGIN_DEPOSIT_01.md` | The static net and K do not select the temporal generator at audited scope. The first-order dissipative form is supplied. |
| `S5_OWNER_RULING_03.md` | A supplied conservative parent derives an underdamped Markov law in its admitted limit, not the inertia-free Level-0 law. No automatic wide-band or overdamped rescue. |
| `L0_ACCESS_BRIDGE_01.md` with its corrections and rulings | The fixed classical net does not generate the later quantum state/access algebra. |
| `L0_LIFT_SELECTION_EVALUATION_01.md`, corrections, ruling 02 | No authorized earned discriminator selects among the surviving physical lifts. Agreement requires a declared per-lift readout and is linear-domain scoped. |
| `CP1_COUPLING_VERDICT_01.md` | Conservation leaves an improvement parameter; dimension/access-dependent restrictions do not discharge geometric coupling. This packet does not reopen that campaign. |
| `GRUT_WORKING_THEORY_SYNTHESIS_01.md` | Net, generator, physical lift, action scale, access and gravitational structure remain supplied at their recorded scopes. |

This pass therefore **does not try to select a quantum algebra by renaming a
classical embedding**. It fixes one ordinary finite quantum control and asks a
narrower question: what interacting information is missing even after that
choice? It is not a rerun of S5 or a proposed new canonical parent.

## 2. Fully specified conventional comparison

Set units of action and reference time explicitly, so that hbar = 1 below.
Restoring units gives `H_physical = hbar H` and energy shift `hbar u`.

The supplied space is `H = C^3 tensor C^3`, with fixed labelled basis
`|n1,n2>`, each ni in {0,1,2}. This is an exact finite three-level inventory,
**not an assertion that finite matrices obey the bosonic CCR**. Define

\[
b=|0\rangle\langle1|+\sqrt2|1\rangle\langle2|,\quad
b_1=b\otimes I,\quad b_2=I\otimes b,\quad
n_i=b_i^\dagger b_i,\quad N=n_1+n_2.
\]

Let Pn be the spectral projector of N for n = 0,...,4 and
P = P0 + P1. The sector dimensions are (1,2,3,2,1). Geometry is a declared
two-site edge. The complete bounded operator algebra, readouts, number
interpretation and low-sector access are supplied.

Use the two-site instance of the record's pinned Laplacian convention:

\[
K=\begin{pmatrix}13/10&-1\\-1&13/10\end{pmatrix},\qquad
\operatorname{spec}K=\{3/10,23/10\}.
\]

This is obtained from `calc/c1_seam.py:build_K(2,0)`, not the literal 23-site
bath used in some canonical runs. Set omega = 5 and u in [0,1]. The autonomous
system Hamiltonian is

\[
H_u=\omega N+\frac u2 N(N-1)
=\omega N+\frac u2\sum_i n_i(n_i-1)+u n_1n_2.
\tag{1}
\]

Its terms are on-site or on the declared edge. All Hu are real, site-exchange
symmetric, number-conserving, positive, have the same unique vacuum ground
state and the same one-particle gap 5. Their operator norm is at most 26.
The interaction is a familiar quartic number interaction, not a generating
principle for GRUT.

Define three loss maps and twelve energy-resolved jumps:

\[
L_1=\sqrt{3/5}\,b_1,\quad L_2=\sqrt{3/5}\,b_2,\quad
L_3=\sqrt2(b_1-b_2),\qquad J_{\mu n}=L_\mu P_n.
\]

For any initial density matrix on this nine-dimensional space, use the supplied
finite GKSL equation

\[
\mathcal L_u(\rho)=-i[H_u,\rho]+\sum_{\mu=1}^3\sum_{n=1}^4
\left(J_{\mu n}\rho J_{\mu n}^\dagger-
\tfrac12\{J_{\mu n}^\dagger J_{\mu n},\rho\}\right).
\tag{2}
\]

All operators are bounded, so the solution exists uniquely for all t >= 0 and
is a CPTP semigroup. No spatial or temporal boundary condition is hidden:
finite open two-site graph, given initial state, forward physical time.
The common rate matrix satisfies

\[
\sum_{\mu,n}J_{\mu n}^\dagger J_{\mu n}
=2\sum_{ij}K_{ij}b_i^\dagger b_j\ge(3/5)N.
\tag{3}
\]

**Every input in this section is POSTULATED as a standard-model control.**
Hu's choice does not generate the sites, quantum algebra, number observable,
clock, vacuum preparation or bath. Their consequences below are DERIVED.
No premise is assigned a frozen L0 information price because no law card is
submitted or evaluated; the specification inventory is explicit instead.

## 3. Exact interaction-completion obstruction

### 3.1 Equality of the entire restricted channel

P is invariant under the forward dynamics: the Hamiltonian preserves number,
and all jumps lower it. On P, `Hu = omega N` for every u and only the n = 1
jumps act. Hence

\[
\mathcal L_u(X)=\mathcal L_0(X)\quad\text{for every }X=PXP.
\tag{4}
\]

Uniqueness of the finite linear ODE gives equality of the channels at **all
times**, not a finite derivative jet. In the basis (vacuum, |10>, |01>), put
`Tt = exp[-(K + i omega I)t]` and write any operator as

\[
X=\begin{pmatrix}x&v^\dagger\\w&R\end{pmatrix}.
\]

The full restricted map is

\[
\Phi_t(X)=\begin{pmatrix}
x+\operatorname{tr}R-\operatorname{tr}(T_tRT_t^\dagger)&v^\dagger T_t^\dagger\\
T_tw&T_tRT_t^\dagger
\end{pmatrix}.
\tag{5}
\]

After the **supplied** carrier-frequency demodulation, its one-particle
amplitude is exactly the Level-0 contraction exp(-Kt).

For any finite adaptive protocol whose preparations and instrument maps send
operators supported on P back into P, every outcome probability is identical
for all u. Proof: insert the common channel (5) between each common instrument;
induct over outcomes and steps. Auxiliary reference systems are also allowed
if the system support stays in P: equality of the linear maps implies equality
after tensoring an identity. Thus increasing lower-sector accuracy, time
resolution or protocol count cannot recover u. A probe creating two particles
is explicitly outside this restricted comparison.

### 3.2 What remains free, rather than being moved by a coordinate change

Fix the vacuum energy and the Hamiltonian's whole P block. For number-conserving
Hermitian completions, the undetected perturbations are exactly

\[
\mathfrak I_P=\{V=V^\dagger:[V,N]=0,\ VP=PV=0\}
=\operatorname{Herm}(3)\oplus\operatorname{Herm}(2)\oplus\mathbb R.
\tag{6}
\]

There are 14 real linear completion directions in this finite inventory.
This is an algebraic count, **not** an information price or a criterion that
selection must reduce dimension. Constants common to the entire Hamiltonian
have already been removed by fixing the vacuum energy.

Even restricting to ordinary quartic interactions leaves nine directions.
Define pair annihilators

\[
A=(b_1^2/\sqrt2,\ b_1b_2,\ b_2^2/\sqrt2),\qquad
V_Q=\sum_{\alpha\beta}Q_{\alpha\beta}A_\alpha^\dagger A_\beta,
\quad Q=Q^\dagger\in\operatorname{Herm}(3).
\tag{7}
\]

In the labelled N = 2 basis (|20>,|11>,|02>), VQ is exactly Q; on P it
vanishes. Thus the map is injective. Q >= 0 makes VQ >= 0 on the full space;
the positive cone has nonempty interior. Real interactions invariant under
site exchange still contain the four-parameter matrices

\[
Q=\begin{pmatrix}a&b&c\\b&e&b\\c&b&a\end{pmatrix}.
\]

These counts concern Hamiltonian completions and their restricted channel
equality. The full energy-conservation construction in section 4 is proved
for the explicit subclass Q = u I, giving (1). It is **not** silently asserted
for every Q. A future theorem fixing arbitrary interacting completions needs
additional information beyond the earned linear sector. Requiring a quasi-free
completion would select u = 0 here by **supplying the absence of interactions**.

It is unnecessary to remove all completion directions to earn a prediction.
Only directions active in the claimed observable matter. Here u is active in
the chosen addition frequency; other observables could be rigid over the same
family. This audit neither restores a unique-universe requirement nor treats
an independently controllable material interaction as something a fundamental
theory must universally hold fixed. It tests only whether the present restricted
data and earned premises already derive that interaction. They do not.

### 3.3 Representation covariance versus different physics

| Transformation | Status | Effect on the claimed observable |
|---|---|---|
| Simultaneously conjugate states, generators, probes and readouts | Passive description change | Measured transition frequencies unchanged |
| Unitarily rotate jump labels within one number-transition sector | GKSL representation freedom | The entire generator is unchanged |
| Apply a number-preserving symmetry to the central interaction in (1) | Physical symmetry if it also preserves the remaining model and fixed interface | N(N-1) is invariant; its coefficient u is not selected |
| Positive rescaling of the clock | Description freedom only if calibration changes with it | chi/kappa is invariant; keeping calibrated K fixed forces scale factor one |
| Change u with the clock, inventory, loss maps and readout fixed | Different autonomous physical generator | The two-particle resonance changes |
| Apply a later pulse or replace the bath/preparation | Different physical history or model | Not a passive equivalence; no uniqueness claim across such changes |

No unitary preserving the declared N can turn Hu into Hu' for u != u':
both are scalar on each number block, so such a unitary fixes each Hu.
More generally, the calibrated energy spectra in section 5 distinguish them.
This is the stopping test the owner requested: the baseline fibre contains
physical completions as well as representation freedom. A claimed relation
must be well defined on the latter and justified across the former; calling
all of them a symmetry would erase real measurable differences.

## 4. Physical admissibility is checked, not inferred from CP alone

For the explicit central-interaction family,

\[
[H_u,J_{\mu n}]=-\nu_nJ_{\mu n},\qquad
\nu_n=\omega+u(n-1)>0.
\]

Therefore

\[
\mathcal L_u^*(H_u)=-\sum_{\mu,n}\nu_nJ_{\mu n}^\dagger J_{\mu n}\le0.
\tag{8}
\]

The system loses energy to a supplied sink; it does not mysteriously conserve
system energy while dissipating. Closed Hu dynamics conserves energy and N.
Both descriptions have a passive vacuum. No finite-temperature KMS/FDT or
Onsager claim is made for the irreversible zero-excitation sink; no such
condition is omitted from a declared finite-temperature equilibrium baseline.

Here is an explicit physical sink, also preventing an abstract-kernel loophole.
For each jump take a fresh two-level ancilla initially in its ground state,
with `HA = nu_n |1><1|`, and define

\[
C_{\mu n}=J_{\mu n}\otimes|1\rangle\langle0|
+J_{\mu n}^\dagger\otimes|0\rangle\langle1|.
\]

It satisfies exactly

\[
[C_{\mu n},H_u+H_A]=0,
\qquad [C_{\mu n},N+|1\rangle\langle1|]=0.
\tag{9}
\]

Each collision conserves total bare energy and total excitation number. Its
support is within the supplied two-site edge plus its attached ancilla. It is
an ordinary bounded finite Hamiltonian, not a claim of relativistic QFT or of
an independently derived lattice. The two-site certificate does not establish
uniform locality of number projectors on longer chains.

For delta > 0 use twelve sequential collisions, each of duration
tau = delta/12, with

\[
W_{\mu n}=\exp[-i\tau(H_u+H_A)-i\sqrt\delta C_{\mu n}].
\tag{10}
\]

Trace out that ancilla after its collision. A round lasts delta. Its resource
ceiling is twelve fresh qubits, ancilla bare gaps <= 8, system energy <= 26,
and interaction norm <= `48/sqrt(delta)` per collision, uniformly in u.
At each fixed delta these are finite and common ceilings, **not identical
microscopic Hamiltonians or free preparations/controllers**. Ancilla gaps in
higher sectors change with u and are explicitly part of the different parent
models. Ground-state ancillas, timing, switching and their replenishment are
supplied. Since the initial ancilla has zero excitation and zero mean coupling,
the mean coupling energy starts at zero and is conserved during that collision;
switching has zero ideal mean work but nonzero control capacity.

For n >= 2, C annihilates P tensor ancilla-vacuum. For n = 1 both gap and
coupling are independent of u. Thus **every fixed finite collision stream has
exactly the same P-restricted channel for all u**, without taking a limit.

The collision stream approximates (2), but is not an exact finite-resource
realization of its exponential time law. There is a whole-domain bound, not
a finite-sample extrapolation: for 0 < delta <= 1 the single-collision channel
is `I + delta D[J] + O(delta^2)`. All odd coupling terms vanish upon tracing
the initial ancilla vacuum. With c = ||J||, its remainder in diamond norm is
at most `(2c)^4 exp(2c) delta^2/24`. Compare to exp(delta D[J]), whose own
quadratic remainder is bounded by `d^2 exp(d) delta^2/2`, d = 2c^2.
For the twelve free/collision factors put

\[
S=2\sup_u\|H_u\|+\sum_{\mu n}2\|J_{\mu n}\|^2\le436.
\]

A conservative splitting bound is `S^2 exp(S) delta^2`. Consequently a single
round differs from exp(delta Lu) by at most C0 delta^2, where the finite,
u-independent C0 is the sum of these disclosed remainder constants.
CPTP contractivity and telescoping give `m C0 delta^2 <= T C0 delta` for
m rounds with m delta <= T. This is a branch-uniform finite-horizon certificate.
The bound is very loose and not a useful device-error estimate. Also the
coupling ceiling diverges as delta tends to zero and the ancilla count grows
as 12T/delta. **No fixed-resource continuum limit is claimed.** The physical
finite-stream obstruction already holds exactly without either limit.

These are standard repeated-interaction constructions. The reviewed primary
collision-model source explicitly discusses the same energy commutator and
the resource cost of a nonzero continuous-time dissipator; see `SOURCES.json`.

## 5. A measurable witness, not a new predicted value

The n-particle energy in (1) is

\[
E_n=\omega n+\tfrac u2 n(n-1).
\]

Using the same calibrated particle-addition probe and phase/frequency reference
in each model, measure the first and second addition frequencies. In units of
angular frequency,

\[
\omega_{01}=5,\quad\omega_{12}=5+u,\qquad
\boxed{\chi=\omega_{12}-\omega_{01}=u.}
\tag{11}
\]

This requires accessible two-excitation transitions and a reference/probe that
can exchange excitations; those resources are declared, not generated. It
does not require breaking total-number conservation of system plus reference.

In the Markov model define the normalized symmetric-mode states

\[
|1_+\rangle=(|10\rangle+|01\rangle)/\sqrt2,\quad
|2_+\rangle=(|20\rangle+\sqrt2|11\rangle+|02\rangle)/2.
\]

For X = |2+><1+|, direct substitution in (2) gives the exact pole

\[
\mathcal L_u(X)=[-3(3/10)-i(5+u)]X.
\tag{12}
\]

The jumps cannot connect its different number blocks because they are
number-resolved; (3) supplies the two anti-commutator eigenvalues. The frequency
shift is therefore present in an operational response, not only an assigned
Hamiltonian coefficient. The first addition pole is `-3/10 - i5`.

The continuum u in [0,1] realizes chi in [0,1] with the *same whole low-sector
channel*. Even after fixing a physical clock, chi/kappa covers [0,10/3].
This proves a non-singleton observable image; numerical fixtures only verify
its implementation.

| u | Same one-particle K and carrier | chi | chi/kappa |
|---:|---|---:|---:|
| 0 | yes | 0 | 0 |
| 1/4 | yes | 1/4 | 5/6 |
| 1/2 | yes | 1/2 | 5/3 |
| 1 | yes | 1 | 10/3 |

A proposed inference of a unique chi from those lower-sector data is refuted
by any pair of these exact models. Conversely, the countermodel itself would
fail if its low-sector channels disagreed or if (11)–(12) failed under its
specified operators. The full matrix-unit and energy checks test those claims.
This is an operationally specified control, not a calibrated experiment or a
novel prediction. Hardware feasibility, preparation errors, pulse bandwidth
and spectral-estimation confidence have not been certified.

## 6. Hostile boundaries and novelty

* **Full dynamics is sufficient.** Once Hu, all jump/collision operators,
  ancilla preparations and interventions are fixed, the probabilities are
  uniquely determined. No cosmic mechanism is needed to choose a pulse or
  material interaction after the model is specified.
* **This is not a faithful completion of the whole nonlinear classical flow.**
  The lower quantum sector is the declared interface. In fact
  `i[Hu,b_i] = -i omega b_i - i u b_i(N-1)`; requiring full-field linear closure
  would detect u and exclude these interactions by an extra assumption.
  Finite capacity also prevents an unrestricted exact CCR/linear-field claim.
* **Statistics/algebra remain supplied.** The finite space admits double
  occupation. A fermionic inventory changes pair sectors. This does not select
  Bose statistics or repair the old lift boundary.
* **Completeness of the measurement record matters.** High-sector probes,
  finite-temperature upward transitions or complete interacting spectroscopy
  can reveal u. Their data are not secretly included in the restricted baseline.
* **Symmetry does not fix an invariant coefficient.** Number conservation,
  reality and exchange symmetry all hold for the whole explicit family. CP,
  positive energy and an energy-conserving completion hold as well.
* **No new mechanism is discovered.** The quartic interaction and its absence
  in the one-particle sector are standard many-body physics. GKSL complete
  positivity and repeated-interaction embeddings are standard. The comparison
  is a transparent null for GRUT's substrate-to-interaction claim, not an
  originality claim for the equations.

## 7. Reproduction and verdict

Run from repository root:

```bash
bash research/interaction_origin01/reproduce.sh
```

The script verifies exact integer/rational counts and energy identities, the
whole restricted channel on all nine matrix units, the full-channel Choi
positivity, total-energy commutators, finite-collision equivalence, physical
resonance shifts, and representation-rotation controls. Inputs, every check,
runtime versions, and convergence diagnostics are in `results.json`; generated
tables are in `TABLES.md`. Hashes are in `MANIFEST.json`. Numerical checks support
the implementation; the all-time/all-u statements follow from the proofs.

**VERDICT: no interaction selector or fundamental GRUT law established.**
The exact scoped obstruction is that a restriction map to the complete
vacuum/one-particle dynamics has a nontrivial interacting completion fibre,
even within positive local, physically realized quantum controls. The current
GRUT record supplies no earned condition that fixes this fibre. This does not
prove that no deeper principle exists, nor suggest that different legitimate
microscopic interactions should have the same response.

The conceptual route closes here: the presently earned lower-sector data cannot
be promoted into a derivation of a physically complete interaction generator.
No further control family is commissioned. One card consumed and killed at
author level, two remain; **Card 2 unopened**. Canonical science, external
evaluator/CR-5 clearance and experimental-lock status are unchanged.

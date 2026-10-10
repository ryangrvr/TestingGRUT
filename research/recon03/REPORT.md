# Bounded local control: bits can move, recoverable records have a horizon

10 October 2026. Author result on the authorized development branch, pending
independent review. **Standard-physics audit; Card 2 unopened.** Scope and
all supplied resources are fixed in [SCOPE.md](SCOPE.md).

## The mathematical answer

In a chain of N memories containing perfect copies of one pointer bit, an
environment-only circuit of at most d layers of disjoint neighbouring gates
cannot remove every record from the predeclared blocks of length 2d+1. Each
such block retains an exactly readable copy. For explicitly defined packing
capacity,

\[
 R_0^{\rm pack}\ge\max\left(1,\left\lfloor{N\over2d+1}\right\rfloor\right).
\]

Moreover, **for even N, the minimum depth permitting R_0^pack=1 is exactly
N/2**. The two endpoint light cones prove necessity; inward nearest-neighbour
CNOT sweeps attain it. This is a circuit result, not a speed limit for every
continuous Hamiltonian. The latter has tails and needs a quantitative
Lieb-Robinson certificate, given conditionally below.

At the same fixed one-layer ceiling, another explicit circuit makes **every
singleton record vanish**, while adjacent pairs retain the bit. Thus locality
protects recoverable capacity in sufficiently large blocks, not each local
bit, a fixed raw detector, or conventional redundancy irrespective of its
fragment convention. The complete reduced S channel is unchanged in all
environment-unitary comparisons.

This answers the owner's finite Wheeler/quantum-Darwinism question in this
declared model. It is not a new law of reality: the restriction follows from
ordinary quantum mechanics and the supplied locality/resource premises.

## 1. Entire domain, preparation and controls

**POSTULATED control inventory.** Hilbert space is C^2_S tensor (C^2)^tensor N_E,
N>=1, with labelled open-chain sites S,E0,...,E(N-1). The pointer is Z_S.
The environment begins in |0^N>. Applying CNOT(S->E0), then
CNOT(E0->E1), ... , CNOT(E(N-2)->E(N-1)), gives the fixed isometry

\[
 W|b\rangle_S=|b\rangle_S|b^N\rangle_E,\quad b=0,1.
\]

That common write costs N sequential neighbour pulses. For an input |+>,
the resulting state is (|0>|0^N>+|1>|1^N>)/sqrt(2). The Born rule,
blank pure memory, labels, distinguished basis, writing circuit and geometry
are inputs. This **constructs copies conditional on those inputs**, not a
derivation of their natural existence or the occurrence of a unique outcome.

Let C_N,d be all E-only circuits U=U_d...U_1 of at most d layers, each a
product of arbitrary two-qubit unitaries on disjoint edges. Identity padding
uses the same ceiling. No fresh ancillas, feedback, postselection or nonlocal
gates. d is a nonnegative integer. Set tau_g=pi*hbar/J for J>0. A gate G has
a Hermitian logarithm A with spectrum in [-pi,pi], G=exp(-iA); choose
H_G=hbar*A/tau_g. Thus every admitted gate has norm <=J and uses tau_g.
Total scrambling time is d*tau_g. Readout has a separately supplied d-layer
budget. Control expressiveness and a calibrated circuit description are
inputs, not free consequences of an arbitrary laboratory's native gate set.

The Hamiltonian in one layer is the sum of its disjoint pair terms. Each
term is bounded by J, and the *global* norm is at most floor(N/2)J. These
are different restrictions. If a different total-norm ceiling J_tot is
independently fixed, pulse time must satisfy
tau_g>=pi*hbar*floor(N/2)/J_tot as well. Every comparator uses the same
ceilings; no parallel control is silently priced at one gate's total norm.

Memory bare Hamiltonians are zero; H_S=hbar*omega*Z_S/2. Both the common
write and E-only controls commute with H_S. Driven memory controls need
not conserve an independently nonzero memory energy. Unitarity supplies
positivity, normalization and ordinary quantum no-signalling. Spatially
local gates supply the specified finite-depth causal cones. No equilibrium
state, bath, KMS/FDT or Onsager claim is made: those conditions do not apply
to this deliberately driven finite register. A relativistic/QFT completion,
controller energy, reset work and native pulse calibration are not derived.

**DERIVED response equality.** Put |e_b>=U|b^N>. Orthogonality is preserved.
On every system matrix unit,

\[
 \mathcal E_U(|b\rangle\langle c|)
 =|b\rangle\langle c|\langle e_c|e_b\rangle
 =\delta_{bc}|b\rangle\langle c|.
\]

This proves identical complete dephasing channels for all U, including
extension by any untouched reference. The visibility is V=0. This audit does
not prove the same perfect-copy bound on every partial-visibility slice of
recon02; perfect initial replication is explicitly load-bearing.

## 2. Fixed operational readout and several distinct record measures

For an interval A, define rho_A|b=Tr_(E\A)|e_b><e_b| and

\[
 D_A={1\over2}\|\rho_{A|0}-\rho_{A|1}\|_1,\qquad
 p_{\rm guess}(A)={1+D_A\over2}.
\]

Equal prior labels b are part of this test. For *each fixed m*, the fragments
are predeclared equal blocks B_j^(m)=[jm,(j+1)m-1], j<floor(N/m). Their
positions, size and number do not change between circuits at a given N,d,m.
The guaranteed chart uses m=2d+1 when it fits. Comparing different m is an
explicit comparison of access classes, not a derivation of an interface.

The binary POVM may depend on the calibrated U but acts only within its fixed
block. A concrete decoder is proved below. In contrast, a *raw Z measurement*
on one qubit has classical distinguishability |rho_0[0,0]-rho_1[0,0]|;
this need not equal D_A.

Define R_delta^pack as the largest number of disjoint nonempty contiguous
intervals whose D_A>=1-2delta, for 0<=delta<1/2. This packing capacity permits
choosing intervals from a fixed declared chart of all intervals. It is not
the averaged quantum-Darwinism statistic N/m_delta. The fixed-block theorem
does not rely on selecting fragments after seeing U.

We separately compute pointer Holevo information

\[
 \chi_Z(A)=S((\rho_{A|0}+\rho_{A|1})/2)
          -[S(\rho_{A|0})+S(\rho_{A|1})]/2,
\]

and actual quantum I(S:A). At D_A=1 the branch supports are orthogonal, so
chi_Z(A)=1 bit and a binary measurement obtains it. Mutual information can
instead reach two bits for access to the whole pure S-E state; it must not
be mistaken for a second classical copy. The finite tables additionally
report R_0.1^avg=N/m_0.1 for the smallest m whose **uniform contiguous-interval
average** chi_Z is at least 0.9. This frozen spatial average differs from an
average over all subsets. No claim about every conventional R_delta definition
is inferred from the packing bound.

## 3. Exact proof of shallow-circuit record capacity

**DERIVED.** For every original site i, set O_i=U Z_i U^dagger. Conjugating
a supported operator by one layer enlarges its support by at most one edge
on either side. By induction,

\[
 \operatorname{supp} O_i\subseteq[i-d,i+d]\cap[0,N-1],\quad
 O_i|e_b\rangle=(-1)^b|e_b\rangle.
\]

Let L=2d+1, take the predeclared blocks B_j=[jL,(j+1)L-1] and centers
i_j=jL+d. O_i_j is supported within B_j and is a Hermitian unitary with
eigenvalues +/-1. Its two spectral projectors form a block-local POVM that
reads b with probability one. Thus D_B_j=1, chi_Z(B_j)=1, and the stated
packing bound follows. If N<L, no such block fits; full-environment
orthogonality still gives R_0^pack>=1. It does not certify a smaller fixed
fragment in that case.

This POVM has an explicit control cost. Start from Z_i, and retain only the
gates intersecting its support as that support is propagated forward layer
by layer. Gates outside this cone cancel in U Z_i U^dagger. The retained
circuit C_i acts entirely in the declared block, has depth at most d and
satisfies O_i=C_i Z_i C_i^dagger. Apply C_i^dagger, then measure Z_i.
That decoder uses at most d local layers, another d*tau_g of available
time, and the known control sequence. It is not a fixed unaided detector.
Decoders on the disjoint certified blocks can be performed together.

**Endpoint strengthening.** O_0 is supported in [0,d], and O_(N-1) in
[N-1-d,N-1], with truncation at the boundaries. When 2d<N-1 these intervals
are disjoint and each contains a perfect copy. Therefore

\[
 R_0^{\rm pack}\ge2\quad\text{if}\quad2d<N-1.
\]

For even N, R_0^pack=1 hence requires d>=N/2. Section 4 supplies a circuit
at precisely that depth. This necessity applies to any admitted circuit,
not only CNOT circuits. No optimality claim for odd N is made.

**Certified error version.** If each actual conditional output differs from
its ideal counterpart by trace distance at most eta_b, contraction under
partial trace and the triangle inequality give
D_A(actual)>=D_A(ideal)-eta_0-eta_1. Applying the ideal decoder yields
p_guess>=1-(eta_0+eta_1)/2. A common unitary approximation with operator-norm
error eta gives eta_b<=eta and D_A>=1-2eta. Readout errors require an
additional independently certified error budget. Simulated roundoff is not
such a physical certificate.

## 4. Hostile ordinary controls under the same operational premises

**Singleton kill in one layer.** In ordered basis 00,01,10,11, let
P=diag(1,i), H=(1/sqrt(2))[[1,1],[1,-1]], and
G_B=CNOT*(HP tensor I). Then

\[
 G_B|00\rangle=|\Phi^+\rangle,\qquad
 G_B|11\rangle=i|\Psi^-\rangle.
\]

Both single-qubit marginals are I/2 for both labels, so all singleton D=0.
Apply G_B simultaneously on (0,1),(2,3),... for even N. Every whole pair
still has D=1. Packing capacity is exactly N/2: at least two sites are
needed for a perfect record and the fixed pairs attain it. Identity at the
same depth-one ceiling instead has N perfect singleton records. This is
an actual state/control change, not a detector relabelling. An arbitrary
two-qubit Hermitian pulse implements G_B within the common norm/time ceiling;
a native-gate implementation would need its own synthesis accounting.

**Raw-detector kill.** Applying H to each memory maps conditional 0/1 to
+/-. Every optimal singleton D remains 1 but every raw Z distinguishability
is zero. Recoverability and an unchanging detector reading are different
physical questions even without correlations.

**Explicit compression.** Divide the chain into blocks of length at most
d+1. Within each block, apply CNOT from the penultimate site into the last,
then move left, one edge per layer; do different blocks in parallel. This
maps |b...b> to |b00...0> in at most d layers. An interval distinguishes
iff it contains a surviving label site, so
R_0^pack=ceil(N/(d+1)). This changes capacity inside the fixed ceiling.

**Sharper balanced compression.** For d>=1 divide into blocks of at most
2d sites and choose a central root. Sweep from both ends inward: on the
left use CNOT(i+1->i) in increasing i; on the right use CNOT(i-1->i) in
decreasing i. Each reset target has a neighbour still equal to b.
Layers on opposite ends are disjoint. For an odd block, serialize the two
last gates sharing the root; its length is at most 2d-1, so depth stays
at most d. Each block leaves b only at its root. Hence

\[
 R_0^{\rm pack}=\left\lceil{N\over2d}\right\rceil,\qquad d\ge1.
\]

Writing R_min(N,d)=min over admitted U of R_0^pack, the universal bound and
the construction quantify the remaining freedom:

\[
 \max\left(1,\lfloor N/(2d+1)\rfloor,
              2\,\mathbf1_{2d<N-1}\right)
 \le R_{\min}(N,d)\le\lceil N/(2d)\rceil\quad(d\ge1).
\]

At d=0, R_min=N. For even N, balanced compression at d=N/2 leaves one
root, attaining the endpoint lower bound exactly. The earlier full
one-sided sweep at d=N-1 is retained as a separately labelled, less efficient
longer-budget control. Neither larger-budget control is a counterexample
to the smaller d baseline. All keep the same system channel, write and N.

**What breaks it immediately.** If an extra fresh blank ancilla at every
site and local SWAP/discard are allowed, all E records can be reset in
parallel while S's channel remains unchanged. The bit moves into the
discarded degrees of freedom. Equivalently the one-site reset channel has
Kraus operators K0=|0><0|, K1=|0><1|. E-only CPTP maps preserve S's reduced
state but need not preserve E's capacity. This comparator changes the
explicit no-ancilla/unitary resource premise; it proves the theorem cannot
be generalized to arbitrary local open-system control for free.

## 5. Continuous-time qualification with a defined error enclosure

Finite local Hamiltonians generally have tails, not exact circuit cones.
Bravyi-Hastings-Verstraete give a primary local-Hamiltonian commutator bound
and a conditional-expectation truncation argument (equation 2). The following
consequence uses an **independently certified** bound for the actual H_E(t):

\[
 \|[U(T)Z_iU(T)^\dagger,W_j]\|
 \le c_0\exp[-(|i-j|-vT)/\xi]
\]

for every one-site unitary W_j, with declared c0>0, xi>0 and v>=0 in
sites/time. Time reversal of the local evolution permits this conjugation
direction too. Their values depend on the interaction bound and graph;
we do not fit them from finite fixtures or claim a device-specific value.

For A=[i-r,i+r] intersected with the chain, define the unital completely
positive conditional expectation

\[
 E_A(O)=2^{-|A^c|}\operatorname{Tr}_{A^c}(O)\otimes I_{A^c}.
\]

It is the product of Haar twirls on all outside sites. A telescoping sum,
norm contraction of each twirl and the certified commutator bound give

\[
 \|O_i-E_A(O_i)\|\le\epsilon_r
 :=\min\left(2,
 {2c_0e^{-(r+1-vT)/\xi}\over1-e^{-1/\xi}}\right).
\]

There are at most two sites at each outside distance, producing the
geometric sum. E_A(O_i) is Hermitian with norm <=1, so
(I +/- E_A(O_i))/2 is a valid local binary POVM. Its branch expectations
are at least 1-epsilon_r and at most -1+epsilon_r. Thus

\[
 D_A\ge1-\epsilon_r,\qquad
 p_{\rm guess}\ge1-\epsilon_r/2.
\]

For delta>0, epsilon_r<=2delta is ensured by
r+1-vT>=xi*log[c0/(delta*(1-exp(-1/xi)))]. Disjoint predeclared blocks
of size 2r+1 then have at least the specified accuracy. This is a conditional
analytic enclosure, not numerically or experimentally certified LR constants.
Unlike the circuit inverse-cone decoder, this POVM's finite laboratory
implementation time is not bounded here. If that readout is not admitted,
the trace-distance result is only a capacity statement. Exact delta=0
does not follow from a nonzero LR tail.

## 6. Reproducible computation and operational test

Run `bash research/recon03/reproduce.sh` at the repository root in the pinned
Python environment in requirements.txt. It produces [results.json](results.json)
and the [generated tables](TABLES.md). Each local unitary and generating
Hamiltonian, seed, interval, fixed-block chart, decoder expectation, channel
matrix-unit error, Holevo value and mutual information is machine-readable.
All intervals are enumerated; dynamic programming computes the largest
disjoint eligible packing. The perfect-record floating threshold is explicitly
diagnostic, not an exact interval certificate. Universal claims rest on the
proofs, not the sampled random circuits.

The controls include identity, Bell encoding, one-sided and balanced CNOT
compression, general seeded bounded Hamiltonians, a fixed-detector failure,
and separately labelled reset/preparation-error boundaries. Results include
all disclosed fixtures, not only successful examples. [MANIFEST.json](MANIFEST.json)
pins this packet's source/results; original recon02 packets and the Card-1
attempt ledger are unchanged.

One finite falsification control uses N=8 after the common write. Predeclare
singleton fragments and blocks [0,2], [3,5], with d=1 and a separate one-layer
decoder budget. Verify the same S dephasing channel using inputs sufficient
for channel tomography; matrix-unit identities are the exact numerical control.
With identity, all singleton D=1. With the stated Bell pair layer, all
singleton D=0, each fixed pair has D=1, and both fixed three-site blocks can
decode the bit with probability one. Random admitted one-layer circuits must
also retain perfect ideal guessing in both predeclared blocks. At a separately
increased d=4 ceiling, the balanced N=8 circuit reduces packing to one.

On a real device, independently certify preparation, locality/crosstalk,
pulse norm and depth, calibration knowledge and decoder/readout errors.
Given branch trace-distance errors eta0,eta1 and a readout outcome error
budget epsilon_read, a certified success below
1-(eta0+eta1)/2-epsilon_read would contradict at least one member of that
**combined model/certificate**, not establish GRUT. No shot count, actual
platform lock or achieved certification is claimed. The superconducting
quantum-Darwinism paper is experimental context with its own supplied gates,
preparations and scrambling settings, not proof that this complete certificate
is currently available to us.

## 7. Premise boundary, novelty and verdict

| Item | Status | Removal or hostile consequence |
|---|---|---|
| Hilbert tensor factors, labelled geometry and pointer | POSTULATED | No identity/geometry/interface generation claimed |
| Pure blanks and common copying interaction | POSTULATED | A product environment uncorrelated with S has no record to protect |
| Initially perfect replicated bit, V=0 | POSTULATED | Partial-visibility or imperfect-copy cases need separate error premises |
| Closed E unitary and no new ancillas | POSTULATED | Local reset/discard removes E records in parallel |
| Fixed d, local terms and pulse ceiling | POSTULATED | Longer admitted depth allows compression to one |
| Calibrated U, block-local decoder and readout budget | POSTULATED | Fixed raw Z can fail even when optimal D=1 |
| Identical S channel under E-only operations | DERIVED | Ordinary partial-trace/unitarity fact, no new response law |
| Fixed-block recoverability and even-N N/2 concentration threshold | DERIVED | Conditional circuit support theorem with constructive controls |
| Continuous-time accuracy enclosure | DERIVED CONDITIONALLY | Requires certified LR bound; no finite native POVM-cost certificate here |
| Naturally generated pointer, copies, identities, Born outcome | NOT ESTABLISHED | No claim or added principle fills these gaps |

The supplied structures are the mathematical assumptions doing the work.
No L0 information-price credit is claimed for this named conventional control.
Premise sensitivity is not a substitute for Appendix-B pricing of a future
authorized candidate. Every quantum/geometry/readout primitive would have to
be priced or genuinely generated in such a candidate.

This is standard quantum circuit theory and information theory, with the
continuous extension a direct application of established local-Hamiltonian
propagation bounds. Predictive-state identification or causal inference cannot
turn a known/calibrated U into an emergent observer. Open-system physics
provides the reset boundary; no Gibbs/KMS/statistical-mechanical novelty is
present. Combining these statements or calling the surviving bit an invariant
does not create an assembled new law.

Wheeler's apparatus-based hypothesis motivates asking about accessible records;
it does not prove that record capacity is the substance of reality. The result
shows an ordinary finite-resource mechanism for **persistence of already
written, decodable records**. It does not show why those records, the pointer,
the partition, the preparation or the access architecture arise.

**VERDICT: retain as a proved conditional standard-physics control, pending
independent review. No new GRUT joint-observable exclusion, generating law,
theory of reality, consciousness theory or experimental lock established.**
Card 1 remains consumed and killed; two cards remain; Card 2 stays unopened.
Canonical state, evaluator/CR-5 clearance and bank flags are untouched.

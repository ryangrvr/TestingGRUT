# Endpoint record allocation and ordinary control closure

10 October 2026. Supplement to recon02, on the same authorized branch. **No
Card 2, generating K or candidate restriction is formulated or evaluated.** This
is a conventional §15.6 NULL-PiT mathematical comparison. Its result is a
conditional obstruction, not new physics, a fourth selector, a universal
impossibility theorem, or a new prerequisite campaign.

## Result, with the condition stated first

If a nonempty family of the recon02 pure conditional memory pairs is closed
under the common environmental unitaries below, with two two-qubit pulse slots
and T/8 of usable time remaining, then at every visibility it realizes, its
endpoint local-record image is the **whole** interval

\[
\{D_F\}=[0,\sqrt{1-V^2}].
\]

Consequently, a relation narrower in D_F cannot remain invariant under this
allowed control freedom. This is a constraint on what an eventual claim would
have to give up or derive; it does not establish that any GRUT family satisfies
the closure condition. It leaves open restrictions on V itself, restrictions
on physical controllability, and relations involving independently measured
costs or specified protocols. None of those is proposed here as a candidate.

## 1. Complete explicit conventional router

Fix v in [0,1], s=√(1−v²). One controlled S-F rotation prepares conditional
memory states e0=|00⟩ and e1=v|00⟩+s|10⟩ from the same |00⟩ preparation.
The S coherence multiplier is v. Let the system-independent F-R unitary be

\[
B_\phi=\begin{pmatrix}
1&0&0&0\\
0&\cos\phi&\sin\phi&0\\
0&-\sin\phi&\cos\phi&0\\
0&0&0&1
\end{pmatrix},\quad 0\le\phi\le\pi/2,
\]

in the ordered basis (00,01,10,11). It fixes 00 and maps 10 to
sinφ|01⟩+cosφ|10⟩. A fully defined bounded Hamiltonian is
H_B=(ℏφ/δt)G, with G=i(|01⟩⟨10|−|10⟩⟨01|), so that Bφ=exp(−iφG).
Its norm is ℏφ/δt. Take δt=T/16, giving norm at most 8πℏ/T, below recon02's
16πℏ/T ceiling. This is one F-R gate, with no fresh memory, conditioning,
postselection or system-dependent feedback. It commutes with H0 on S and with
the excitation number of F+R; both memory bare Hamiltonians remain zero, as
in the original P. No thermal bath, KMS/FDT or relativistic completion is added.

After Bφ, the conditional F states are |0⟩⟨0| and

\[
\rho_{F|1}=\begin{pmatrix}
v^2+s^2\sin^2\phi&vs\cos\phi\\
vs\cos\phi&s^2\cos^2\phi
\end{pmatrix}.
\]

Their difference is traceless, with eigenvalues ±D_F. Thus

\[
D_F^2=s^2x(v^2+s^2x),\qquad x=\cos^2\phi.
\]

At x=0, D_F=0; at x=1, D_F=s. For s>0, the polynomial is continuous and
strictly increasing on [0,1], except that its derivative is zero at the single
endpoint x=0 when v=0. It still takes every value in [0,s²]. For any desired
d in [0,s], the exact choice is

\[
x=\frac{\sqrt{v^4+4d^2}-v^2}{2s^2}.
\]

For v=1, the only allowed d is zero and the router is immaterial. The operation
preserves the overlap exactly, not merely its modulus. Partial trace over FR
therefore preserves the full reduced S channel on every input, including inputs
entangled with an external reference. This is a physical operation on the
memories with a fixed F readout, not a detector relabelling or an R1 quotient.

At v=3/5, φ=0 gives D_F=4/5 and φ=π/2 gives D_F=0. The intermediate value
D_F=3/5 is obtained with

\[
x=(\sqrt{981}-9)/32.
\]

The same initial controlled rotation suffices; only the subsequent ordinary
memory router changes. Unlike recon02's product examples, this control does
not require retuning the original S-memory coupling to span the slice.

## 2. Arbitrary branch pairs: exact transitivity

For any normalized pure pair e0,e1 in C⁴, set q=⟨e1|e0⟩ and v=|q|. If v<1,
the two vectors

\[
f_0=e_0,\qquad f_1=(e_1-q^*e_0)/\sqrt{1-v^2}
\]

are orthonormal. Complete them to an orthonormal basis f0,...,f3. Map this
basis unitarily to (|00⟩,|10⟩,|01⟩,|11⟩). The resulting common unitary Wc
maps the pair to |00⟩ and q*|00⟩+√(1−v²)|10⟩. For v=1, e1=q*e0, and a
one-dimensional isometry extended to a unitary gives the same conclusion.
No branch selection or intervention-dependent map is used: one fixed Wc acts
on both alternatives. The original complex q, hence the whole S channel, is
preserved. The formula in §1 uses |q|² for a complex overlap and is unchanged.

Every finite-dimensional Wc has a Hermitian logarithm A with spectrum in
[-π,π], Wc=exp(−iA). Implement it on neighbouring F-R in one pulse of length
T/16, with norm ≤16πℏ/T. Follow it by the one routing pulse. This uses two
two-qubit pulse slots and T/8, not extra ancillas or a larger memory. The
existence of this logarithm is a control assumption of the declared P, not a
claim that an actual device synthesizes every such pulse with this cost.

## 3. The conditional family theorem

Let C be a family of the original pure-pair protocols. Suppose a member at a
given q has at most four occupied two-qubit pulse slots and enough spare time,
and C also admits its Wc and Bφ extensions under the same ceiling. These two
closure requirements are explicit logical hypotheses. By §2 it contains a
canonicalized pair; by §1 its observable image at the same q contains every
D_F in [0,√(1−|q|²)]. Recon02's upper bound gives the reverse inclusion.
Therefore that slice is exactly the full interval.

If these assumptions hold at every realized visibility, then

\[
O(C)=\bigcup_{v\in V(C)}\{v\}\times[0,\sqrt{1-v^2}].
\]

If V(C)=[0,1], O(C)=S_P. In this case no new joint restriction on (V,D_F)
exists. If V(C) is smaller, the restriction is on visibility; the record
allocation remains unrestricted at each surviving visibility.

This does **not** make all resource-limited protocols closed under arbitrary
extra gates: a six-slot protocol may have no reserve, and a long protocol may
have no time reserve. The theorem asserts closure only where it is assumed
and operationally fits. The canonical seed plus router uses just two slots,
so its full interval already fits the original ceiling directly. In the
arbitrary-pair theorem, preparation/analysis pulses still count against the
eight one-qubit slots. The two added F-R operations add no one-qubit slots.

The original report used “entangling two-qubit gates” while already counting
SWAP in its synthesis. Here “two-qubit pulse slots” makes that intended
operation count explicit; SWAP itself is not an entanglement-generating gate.
The budget is not increased. If a device instead defines a slot as one native
CNOT or another restricted gate, the arbitrary-Wc one-slot implementation is
not certified by this packet: it requires its actual synthesis and resource
cost. The explicit router has its stated G Hamiltonian; whether G is available
is also an independently testable physical control premise.

## 4. Approximate operations: quantitative measurement boundary

If the implemented common environmental unitary is Wtilde with
‖Wtilde−W‖op≤η, both normalized conditional output vectors differ by at most
η in norm. Each corresponding pure-state trace distance is at most η.
Contractivity and the triangle inequality then give

\[
|\widetilde D_F-D_F|\le2\eta.
\]

This bound is independent of dimension. V and the S channel remain exactly
unchanged under any genuinely environment-only unitary, including an imperfect
one. System leakage, residual S-memory coupling and nonunitary control errors
are separate applicability errors, not absorbed in this claim. If two compared
implementations each have the bound η, their nominal D_F separation can shrink
by at most 4η. At v=0.6 the routed endpoints remain separated whenever η<0.2;
tomography uncertainty must also be subtracted. This is a conventional error
certificate on specified operations, not an experimentally measured η.

## 5. What this settles, and what it does not

Known dilation freedom under environmental unitaries is already part of
standard quantum channel theory [1]. The elementary branch-pair proof here
fixes the particular finite problem and finite-control qualification, rather
than promoting that theory to GRUT. Published photonic quantum-Darwinism work
also investigates environmental correlations and their effects on records [2];
it does not certify this router or this pulse ceiling.

A narrower local-record relation can be meaningful only if its family fails
one of the stated closure assumptions, or if the relation is restricted to a
specified stage/protocol rather than all allowed endpoints. Supplying a routing
prohibition would be a premise. Deriving a physically measurable prohibition
would be a different scientific result; none has been derived here. Changing
F to “whatever fragment contains the information” changes the observable and
does not beat recon02's fixed-interface comparison. Declaring an environment-
only operation part of the readout equivalence class would erase the very
physical record allocation being tested; this packet grants no such interface.

The result does not forbid a theory of reality, a generated observer, a
control-sensitive relation or physics beyond ordinary unitary controllability.
It rejects an unconditional endpoint allocation claim **under the stated
closure assumptions**. It does not justify saying a working GRUT theory is
almost established. No selector, generating law or experiment is certified.

## 6. Reproduction and disposition

Run `bash research/recon02/reproduce_routing.sh` from repository root. The script
checks G's exponential, the whole S matrix-unit channel, exact rational
polynomial controls, inversion fixtures, arbitrary-pair canonicalization,
complex overlaps and degenerate endpoints, and calibrated unitary-error
fixtures. Inputs, results and content hashes are machine-readable in
[routing_results.json](routing_results.json) and
[ROUTING_MANIFEST.json](ROUTING_MANIFEST.json). Floating-point diagnostics
support the proofs; they are not interval-certified numerical calculations.

**VERDICT: CONDITIONAL CONTROL-CLOSURE OBSTRUCTION; NO GRUT LAW.** Author proof
and conventional synthetic controls, pending actual independent review. Card 2
unopened; one killed frozen slot consumed, two remain. Evaluator/CR-5 external
clearance remains unmet. No canonical file, bank flag or Card-1 ledger is changed.

## Primary sources

[1] D. Kretschmann, D. Schlingemann, R. F. Werner, *The Information-Disturbance
Tradeoff and the Continuity of Stinespring's Representation* (2006),
[primary manuscript](https://arxiv.org/pdf/quant-ph/0605009), introduction and
§III.A-B: environmental unitary freedom and minimal/nonminimal qualifications.

[2] M. A. Ciampini et al., PRA 98, 020101(R) (2018),
[primary paper](https://pureadmin.qub.ac.uk/ws/files/161849361/Experimental_signature_of_quantum_Darwinism_in_photonic_cluster_states.pdf),
engineered environmental correlations. No original routing theorem is claimed.

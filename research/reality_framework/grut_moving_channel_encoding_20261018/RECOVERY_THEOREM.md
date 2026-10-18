# Recoverability with the declared color output

Tier: author-derived, independent mathematical review pending. The exact
B2-IS decomposition and the triplet observation in B2 diagnosis are inputs,
not new results. This note optimizes recovery and applies that structure
to a separately named protocol; it proves no delivery lower bound.

## Output access and what the environment can record

Write P_+,P_- for the rank3/rank1 SWAP projectors. For the actual law,
the successful dilation is V=P_+ tensor chi_+ +P_- tensor chi_-.
Put u=||chi_+||^2,v=||chi_-||^2,c=<chi_-,chi_+>. These vectors include
the actual receiver purification if its initial state is mixed. Tracing
the environment gives the already derived B2-IS operation E.

For pure input with reference and triplet weight x, the unnormalized
complementary output is x|chi_+><chi_+|+(1-x)|chi_-><chi_-|.
Off-diagonal parity coherence disappears from this complement because
the output color sectors are orthogonal. The environment can record
parity; it need not record an individually identifiable particle label.
Input-dependent success u!=v also lets the flag filter that parity.

The purified conditional parity-vector distinguishability is
sqrt(1-|c|^2/(uv)), when u,v>0. For a MIXED physical receiver, tracing
its inaccessible purifier can reduce physical-environment distinguishability.
Do not equate this purified formula to an actually accessible record.
Even orthogonal final states need not carry an accessible classical label.

Identity and SWAP histories reaching the SAME final environment interfere.
Their alternatives are not measured histories. The current decoder traces
spatial, tail, fuel and receiver registers after extracting the two colors;
it has no flag saying which history occurred. Reading those registers,
feed-forward or an inverse filter would change the decoding protocol.

## Optimal deterministic color-only recovery: a quantitative limit

Allow, as a GENEROUS upper envelope, any deterministic CPTP operation D
on the two output colors, with fresh ancillas but no environment access.
This is broader than the original fixed left/right decoder. No extra
postselection is included. For a four-dimensional maximally entangled
input with a reference, success is q=(3u+v)/4. If q>0, the maximum
conditional entanglement fidelity to that input is exactly

    F_max=(9u+v+6|c|)/(4(3u+v)).                         (1)

Proof: write Kraus M_k=a_k P_+ +b_k P_- for E, with sums |a|^2=u,
|b|^2=v,a conj(b)=c. For recovery Kraus D_l put X_l=Tr(P_+D_lP_+),
Y_l=<singlet|D_l|singlet>. The unnormalized Choi fidelity is

    [u sum|X_l|^2+v sum|Y_l|^2
                         +2Re(c sum X_l conj(Y_l))]/16.

CPTP completeness and Hilbert-Schmidt Cauchy-Schwarz give sums <=9,<=1
and absolute cross sum <=3. The bound is attained by
D=P_++exp(i arg(c))P_- (any phase if c=0). Divide by q only now.
The same proof permits arbitrary ancillas in D. Its scope is this Choi
test, not an optimization of every possible communication protocol.

Perfect deterministic recovery of ALL four input dimensions, with a
state-independent nonzero successful channel, is possible iff

    u=v=|c|>0.                                           (2)

Indeed 1-F_max has numerator 3(u+v-2|c|), which is nonnegative and
vanishes only under (2), since |c|<=sqrt(uv). The phase unitary then
corrects E. The UNCHANGED decoder needs c=u=v>0, not merely (2).

A necessary condition for 1% trace-distance accuracy for every input
and reference, even with that generous recovery, is

    |c| >= (36u+37v)/75.                                (3)

This follows from F_max>=.99 (testing the maximally entangled projector).
It is necessary, NOT sufficient for the worst-input task.

A second necessary condition is |c|>=49(u+v)/100. The two orthogonal
equal-weight triplet/singlet superpositions have conditional output
distinguishability 2|c|/(u+v). Contraction and the triangle inequality
force worst-input recovery error >=1/2-|c|/(u+v). This can be stronger
than the Choi test. Neither condition bounds the actual moving entries.

In the balanced case u=v=p>0, the optimal WORST-INPUT conditional trace
distance is exactly (1-|c|/p)/2. Phase correction attains it. For a lower
bound choose two orthogonal superpositions of one triplet and the singlet;
their output distinguishability is |c|/p. CPTP contraction and the triangle
inequality force either recovery error to be at least (1-|c|/p)/2.
Thus the 1% balanced criterion is |c|>=.98p; the Choi test alone is weaker.

## MC-E1-T: an explicit lower-input-dimension protocol

Use the SAME Hamiltonian, spatial preparation, fuel loading locations,
successful event, trace-out and time windows as the moving experiment.
Change only the supplied message encoding to

    W|0>=|00>, W|1>=(|01>+|10>)/sqrt(2), W|2>=|11>.

This is one joint qutrit in the two fuels, not arbitrary two-qubit data.
It requires preparation of those joint color states, already a restriction
of the arbitrary-AB input resource. The logical output remains distributed
in the same two color ports. W is the observable dictionary, not access to
discarded geometry. No separate local receiver is claimed to read it alone.

Because SW= W, complete evolution factorizes as U_+xi tensor W|message>.
For every time, actual conditional receiver and external reference,

    E_R(W rho W*) = u W rho W*,
    failure reference = (1-u)rho_Z, retention probability=u_T.  (4)

Thus the complete logical instrument is exactly qutrit identity on success
and erasure otherwise. All geometry excursions, conversions and coherent
returns are present in u,u_T. There is ZERO conditional logical error when
u>0, without a small crossing estimate. This is inherited symmetry used as
an explicit protocol, not a new protection mechanism or operating pass.

The largest subspace noiseless for EVERY member of the full nontrivial
interaction algebra span{I,SWAP} has dimension3. A universally scalar SWAP
action must lie in an eigenspace (ranks3 and1). Particular actual channels
can transmit all four dimensions, so this is not their universal capacity
ceiling. Three logical levels also is not a Shannon-capacity calculation.

The operating obligations now reduce to actual u>=.1 and u_T>=.25 over
their original windows. Neither bound is established. If u=0, no message
arrives. A protected alphabet cannot supply matter retention or output.
Original arbitrary-two-qubit operation remains a separate unresolved task.

## Relation to established methods

This uses standard interaction-algebra/noiseless-subspace reasoning, not
a new quantum principle: Knill, Laflamme, Viola, PRL84,2525 (2000),
https://arxiv.org/abs/quant-ph/9908066; Knill and Laflamme, PRA55,900
(1997), https://arxiv.org/abs/quant-ph/9604034. The model-specific work is
the pinned algebra, output access and finite diagnostic below. No agreement
with these mathematical methods is empirical evidence for GRUT.

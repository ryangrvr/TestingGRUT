# B2-SO: joint scaling error controlled; source transfer unresolved

**Decision:** the unchanged fixed-path law now has an author-derived
two-source comparison whose error stays below 0.001 throughout the frozen
growing use interval for L in 24N, L>=2016. The actual outgoing source
probability and late fuel-return bound remain unproved. The operating gate
is therefore **mathematically unresolved**, not passed or refuted.

## What changed

The new argument connects the original double-fuel seed to the SAME joint
receiver projector. It does not just strengthen a bulk band description.

Let u(t) be the exact half-line one-source evolution of Pi f, without
normalizing or physically preparing that filtered state. Cut its two
copies before the center and embed them into the full two-resource Hilbert
space. The physical Hamiltonian retains every constituent position, fuel,
field, interaction, reverse conversion and return. The mathematical
reference differs from the original evolution by

    delta_L(t)<=12/27613+1500 exp(-L/105)<0.001,
                    0<=t<=L/(6J), L>=2016.

The proof constructs an exact distance-weighted invariant graph. All
virtual charged states remain; no finite-order phase approximation is
propagated for O(L/J) time. High-band input norm is paid once. The source
packets' boundary-shell drive is then integrated in exact Duhamel, including
the effect of every later physical excursion and coherent return.

Define a_L as u's co-pair probability in the frozen left receiver, and
m_M as its matter probability inside the half-host cut. The ORIGINAL
physical probabilities obey

    p_out^(L)(t)>=[a_L(t)-delta_L(t)]_+^2,
    p_q=2^(L)(t)>=[m_M(t)-delta_L(t)]_+^2.

Thus multiplying one-source successes is justified only inside this
controlled reference, with its complete error paid in the joint event.
The same physical outgoing projector has the exact endpoint two-qubit
identity channel on success, including entangled inputs. All failures
remain erasures; information quality does not fill the probability gap.

## What is still missing

SOURCE_ESTIMATE.md expresses the remainder with the actual UNFILTERED
one-source fuel state v(t)=exp(-ith)f:

* Its co-pair probability in the prescribed receiver must have a uniform
  lower bound over I_L=[L/(8J),L/(6J)].
* Its fuel return amplitude must stay sufficiently small from 2/J through
  L/(6J), not only near the initial production peak.

For example, receiver probability >=0.34 and fuel-amplitude modulus <=2/3
would suffice after reconstruction and comparison errors. Neither inequality
has been established. The exact source Green function is known, but its
spectral measure and outgoing spatial amplitudes still need quantitative
estimates, including trapped and slow components. No source probabilities
were numerically evaluated this cycle.

There is a narrower original-preparation corollary: combining the comparison
with B1's existing launch bound gives q=2 probability >0.323761 for
2<=Jt<=4 at these lengths. This counts two opposite-charge pairs' worth of
matter and allows separated constituents. It is not co-pair arrival in the
distant receivers or the full growing retention window. It inherits B1's
pending proof review.

## Verification and scope

The frozen rational inequalities passed. Matrix controls at L=4,6,8 checked
the embedding, inherited cut energies, boundary-only defect and derivative
identity; the data replay was byte-identical. Combined workers used 0.641
seconds and peak RSS 44.54 MiB. No trajectory, new candidate, source layout,
parameter change or favorable replacement endpoint was tested.

All-size proofs are author-derived pending independent review. Controls
share B1 code and use ordinary double matrix arithmetic; they are not
outward-rounded certification or external validation. B2 remains failed
at its original eight-event endpoint. A1 and the old solver branch remain
closed. This remains a dilute, nominal-parameter model on prepared fixed
geometry; no moving-geometry, scattering-operation, extensive-world or
distinctive empirical claim follows. The bounded cycle stops at the
explicit source estimate.

Proof: PROOF.md. Remaining estimate: SOURCE_ESTIMATE.md. Review and
challenge: REVIEW_PACKET.md and CHALLENGE.md.

# Actual right-boundary action: an exact compatibility theorem

**Tier:** author-derived for the pinned unchanged B1 law; independent
mathematical review pending. No full-window operating result is asserted.

## 1. Sector, preparation and quantifiers

Fix a supplied recorded simple word w=(w_0,...,w_(m-1)) among N labels,
3<=m<=N, with q=a_w. Put h=w_(m-1), p=w_(m-2), r=N-m and b=a_contact,
c=b* b. The receiver is the ACTUAL 64N^2-mode open chain, with total boson
number m-1, not a smaller substitute. Let sigma_w be its actual conditional
state after recording. No value of Tr(sigma_w c) is assumed.

There are two excited fuel slots, at w_0 and h, initially vacuum matter
and electric fields. Their arbitrary two-qubit/reference state is supplied
independently of sigma_w. Both directions of conversion remain on. The
conversion program stays on these event labels if the head changes.

Work in the EXACT invariant sector selected by that preparation: total
F+N_+=F+N_-=2 and every NONPROGRAM fuel slot is vacuum. B1 has fuel
registers at all events, but only the two programmed labels can convert;
fuel has no hopping term. Thus inactive fuel vacuum is preserved by the
unchanged law. No reachable state is removed by this sector restriction.

Let P fix this entire recorded word and q=a_w, retain EVERY physical
two-resource matter/fuel configuration on it and EVERY receiver state with
m-1 bosons. Thus F+N_+=F+N_-=2, capacity-one fields obey Gauss, and free
events have no matter or excited fuel. P does not require that fuels remain
excited. It contains the whole fixed-path B2 experiment, not just its seed.
Within that invariant sector write Q=1-P, H=H_full, A=P H P, B=Q H P,
C=Q H Q. These are mathematical
compressions; no physical transition is removed.

On P, the exact compression is

    A = H_m tensor I_receiver + I_matter tensor dGamma(h_receiver)
                         + sum_i delta_i(a_i-R_i) I.             (1)

H_m is precisely THEOREM.md in the parent consolidation, including the
endpoint degree correction. All geometry-changing entries leave P. The
preparation is psi_0=|two fuels, w,q=a_w> tensor eta_w, where eta_w purifies
sigma_w and may include the input/reference purification.

For comparison with B2-RA, m=L must already belong to its frozen length
family. This note neither selects a new record nor assigns a probability
to obtaining that length. The identities apply to any given record and
retain the actual conditional receiver, including mixed states.

## 2. Enumerate all geometry-changing outputs

Let n_x=n_+,x+n_-,x, epsilon=1/4, and let G_h project onto matter/fuel
vacuum at h and zero electric field on (p,h). Distinguish these four maps:

| Primitive | Required gates on P | Coefficient/operator | Target |
|---|---|---|---|
| Assisted append of each free v | New child has empty matter/fuel, q_v=0, new field zero; no old-head emptiness | -(g/N)b* | wv, q_v=1, one more receiver boson |
| Bare append of each free v | Same empty-child/field gates; no storage raise | -(epsilon/sqrt(N)) n_h | wv, q_v=0, receiver unchanged |
| Assisted pop of h | G_h, q_h=1, contact annihilation nonzero | -(g/N)b G_h | w^-, q_h=0, one fewer receiver boson |
| Bare pop of h | G_h | -(epsilon/sqrt(N)) n_p G_h | w^-, q_h=1, receiver unchanged |

Each retains its exact adjoint. A bare append from P creates an allowed
storage-error state; a bare pop leaves a stored bit on a free event. These
are physical Q states, not discarded defects. The receiver constraint
Q_storage-1=n_receiver remains satisfied in every case.

Different free children have different words. Assisted and bare outputs
have different storage/receiver sectors; append and pop have different
lengths. Consequently their images are mutually orthogonal on P. No matter,
fuel, field or receiver-only primitive changes the word or leaves P.
The exact positive operator identity is therefore

    B*B = (g^2/N^2)[r(c+1)+c G_h]
                  +(epsilon^2/N)[r n_h^2+G_h n_p^2].             (2)

Equation (2) is an operator equality, including correlations in arbitrary
P states. It is not an incoherent replacement of the time evolution.
Orthogonality is used only for these single-step outputs. In particular,

    B*B >= (g^2 r/N^2) P.                                      (3)

For m<N, B is injective. NO nonzero invariant subspace of H can be
contained wholly in this fixed-word endpoint sector: invariance would
require B psi=0. This rules out an EXACT embedding there, not an approximate
or dressed embedding extending into Q.

## 3. Action on the actual two-fuel seed

Initially n_h=0 and G_h=0. Pop is blocked by the excited head fuel; bare
motion counts matter rather than fuel and vanishes. Thus

    B psi_0 = -(g/N) sum_(v free)
                   |wv,q+e_v; same two fuels,E=0> tensor b* eta_w,  (4)
    ||B psi_0||^2 = (g^2 r/N^2)[1+Tr(sigma_w c)] =: gamma_w^2.     (5)

Colors/reference are acted on as identity. There is no cancellation between
the orthogonal child words and no favorable receiver replacement. Since
0<=Tr(sigma_w c)<=m-1,

    g sqrt(r)/N <= gamma_w <= g sqrt(mr)/N <= g/2.                (6)

For finite N the exact leakage probability has the short-time expansion

    Tr[Q exp(-itH) rho_0 exp(itH)] = gamma_w^2 t^2 + O(t^3).     (7)

The remainder is NOT bounded through the B2 operating interval here. An
unknown actual contact population cannot be replaced by its maximum to
claim observed order-one leakage, nor by zero to claim protection.

The initial growth is storage-consistent. On matter vacuum its diagonal
detuning is delta_v plus the actual receiver creation dynamics:

    [D_receiver,b*]=a*(h_receiver e_contact).

There is no 64-unit consistency penalty for this assisted output. The
large penalty of a bare storage defect therefore cannot block (4).
From a consistent record, bare append has storage cost 48 and bare pop
cost 80 in the actual relative-energy convention. On the co-pair witnesses
below these become 48+A_degree and 80-A_degree after degree changes.
They are energy costs, not exact prohibitions or a proven spectral gap
for the full coupled boundary problem.

## 4. Conversion does not lock the boundary

The on-event term -J exchanges right fuel and the co-pair at h, retaining
its adjoint. It stays inside P until geometry changes.

| Right site state | Assisted append | Bare append per child | Pop |
|---|---|---|---|
| Excited fuel, no matter | Active | Zero | Blocked |
| Co-pair, empty fuel | Active | -2epsilon/sqrt(N) | Blocked |
| One constituent, empty fuel | Active | -epsilon/sqrt(N) | Blocked |
| Empty matter/fuel | Active | Zero | Assisted/bare reverse maps can act |

The last row retains the contact and parent-occupation gates in section 2.
After conversion the bare growth norm on an old-head co-pair is exactly
2epsilon sqrt(r/N), independent of the receiver. This conditional fact
does not assign a population to that co-pair.

After any append the old fuel label h becomes interior. Let I_v append an
empty event to a fixed matter state. On the embedded old state,

    (H_matter^(wv) I_v-I_v H_matter^w) phi
       = A_degree d_h I_v phi
              -jp T_+(h->v) phi -jm T_-(h->v) phi,              (8)

where A_degree=13/46080, jp=1/64, jm=1/96, d_h is co-pair occupation,
and T includes the required field shift. On a co-pair at h these three
outputs are orthogonal, with squared norm A_degree^2+jp^2+jm^2. They have
no distant-head suppression once this geometry exists. Mass, density on
the new empty child, fuel energy and static conversion add no other drive.
Thus the grown geometry is NOT a decoupled spectator of the old source.

These changes are reachable from the actual seed. For a specified free v,
the minimal target with plus at v, minus at h and right fuel empty has:

    assisted q_v=1 target:  H^3 coefficient
                      -2g J jp/N times b* eta_w;               (9)
    bare q_v=0 target:      H^3 coefficient
                      -2epsilon J jp/sqrt(N) times eta_w.       (10)

In (9) growth/conversion have two allowed orders, then outward hop; in
(10) conversion must precede bare growth. The minus-outward versions
replace jp by jm. Both targets obey Gauss and capacity. Minimal-word
identities are not operating-time transition probabilities.

## 5. Fully assembled record is not an exact exception

If m=N, initial assisted growth is absent and B psi_0=0. This does not
establish invariant evolution. Convert the right fuel; move both
constituents to p (two allowed hop orders); bare-pop the now-empty h.
The target is w^- with unchanged q=a_w, the left fuel still excited and
the right co-pair at p. The actual receiver is UNCHANGED. Its minimal
four-step coefficient is

    2 J jp jm (2epsilon/sqrt(N))
                   = J jp jm/sqrt(N),                          (11)

tensor identity on eta_w. No contact occupation is needed. The storage
error is allowed even though costly. An assisted-pop counterpart has
2g J jp jm/N times b eta_w and may vanish on contact vacuum. Hence the
entire two-fuel trajectory is not confined exactly to P even at m=N.
This is not a claim of appreciable full-assembly preparation or late loss.

## 6. Exact comparison route and its limitation

Let H_0=A direct-sum C be the mathematical block reference. Its initial
evolution stays in P and realizes H_m with the genuine receiver. Duhamel
gives the valid complete-state estimate

    ||(exp(-itH)-exp(-itH_0))psi_0||
          <= D_abs(t):=integral_0^t ||B exp(-isA)psi_0|| ds.     (12)

The integrand squared is exactly (2) in the reference state. In this
reference the actual receiver evolves with h_receiver and factors from
matter, so its contact expectation and the endpoint occupation functions
could be calculated without changing preparation. However (3) gives

    D_abs(t)>=t g sqrt(r)/N.                                   (13)

For an accepted m>=N/2, r>=1 and the B2 upper time T=m/(6J)=640m,
the right side is >=320g sqrt(r)>=288 sqrt(r). Therefore this ABSOLUTE
drive certificate cannot establish a 1/1000 state comparison over the
full window, even with perfect knowledge of the reference populations.
This is a limitation of that positive majorant, NOT a lower bound on the
actual error: coherent cancellation, energy denominators and returns can
make the actual difference much smaller. Other, phase-sensitive estimates
remain possible. At m=N the r baseline vanishes and (11), not (13), is the
exact-invariance obstruction.

## Conclusion

Exact Route A fails for the supplied endpoint sector: immediately for
m<N and through a minimal retreat word for m=N. Approximate operation is
not decided. CHANNEL_PROBLEM.md gives the complete retained-memory and
output reconstruction needed for Route B, without treating survival in P
as the moving experiment's success event. B2-RA's fixed-path result is
unchanged; no successor or frozen physical threshold is introduced.

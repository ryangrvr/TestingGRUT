# The changed boundary-source problem, with no deleted environment

**Tier:** exact author-derived reduction; quantitative operating bounds
unresolved. Independent review pending. Definitions A,B,C,P,Q are those
of BOUNDARY_PROOF.md, not a one-particle surrogate.

## 1. Exact history-dependent dynamics

For the original preparation psi_0 in P, write psi_P=P exp(-itH)psi_0.
Solving the Q equation and substituting gives

    psi_Q(t)=-i integral_0^t exp[-iC(t-s)] B psi_P(s) ds,
    i d psi_P/dt = A psi_P(t)
          -i integral_0^t B* exp[-iC(t-s)]B psi_P(s) ds.          (1)

The initial memory value is exactly B*B in BOUNDARY_PROOF equation (2).
All subsequent growth, retreat, storage errors, receiver propagation,
reverse conversion, charged virtual states and coherent returns reside
in C and the memory kernel. There is no Markov or independent-emission
assumption. The fuel program remains fixed on its recorded labels.

A equals the fixed-path matter Hamiltonian plus actual receiver and a
constant, but the kernel is NOT scalar or a receiver-only operator.
Equation (8) of the boundary proof identifies the new diagonal and
outward matter action once the old source becomes interior. No gap for
this whole complement or small real-axis resolvent is supplied.

## 2. Use the same moving readout as the parent packet

Where m=L has B2-RA's frozen length and time ranges, use the same X_L,Y_L
and original recorded labels. Success requires the original prefix,
bonds and its storage bits intact, both selected fuels empty and co-pairs
at one event in X_L and one in Y_L. Decode their plus colors in recorded
left/right order, trace every tail/receiver register, and flag every
other outcome as erasure. No intermediate measurement is introduced.

Growth beyond the original head is allowed by this event. Thus R P and
R Q can BOTH contribute. The P memory equation alone does not equal the
complete output probability: psi_Q must also be reconstructed in (1).
Because R is a configuration projection commuting with P, the final P/Q
parts are orthogonal; interference WITHIN each part remains physical.
Retention uses the same joint prefix/storage survival and k=2 event
specified in the parent moving-composition packet.

No new geometry selection or probability of obtaining an appropriate
record is established. This remains a conditional experiment after the
supplied record and fuel loading, not autonomous preparation.

## 3. Actual successful information instrument

With exactly two color carriers, choose an encoding by their occupied
register order in the current word. At the root place its fuel before
its plus register; at the other selected event place fuel after plus.
Every color-blind primitive acts as identity or SWAP S on the two input
colors. Geometry and receiver registers remain in the spatial factor.
Consequently H=A_c tensor I+B_c tensor S for Hermitian spatial operators;
this A_c is NOT the compression A above. Set

    H_+=A_c+B_c, H_-=A_c-B_c, U_+=exp(-itH_+), U_-=exp(-itH_-).

Purify the actual conditional receiver in the initial spatial state xi_w.
Let chi_+=R U_+ xi_w, chi_-=R U_- xi_w and Pi_+,Pi_- be the color
triplet/singlet projectors. The EXACT complete successful instrument is

    E_R(rho)=sum_(s,s' in {+,-}) <chi_s',chi_s> Pi_s rho Pi_s',  (2)

also tensor identity on any external reference. It follows that

    p_R(rho)=p_+ Tr(Pi_+ rho)+p_- Tr(Pi_- rho),
    p_+=||chi_+||^2, p_-=||chi_-||^2,
    coherence coefficient c_R=<chi_-,chi_+>.                    (3)

Equivalently a=(chi_++chi_-)/2, b=(chi_+-chi_-)/2 give the identity/SWAP
Gram instrument archived in B2. Exact identity on successful arbitrary
inputs requires chi_+=chi_- (equivalently b=0), not merely large p_+ and
p_-. A collective triplet code would be a different task.

The interior-fuel ten-operation exchange word is still allowed. For
three consecutive events a,h,v with two co-pairs at a and v and empty
fuel at h: move both a->h, absorb; move both v->h then h->a; emit at h,
then move both h->v. Every link remains capacity one; all reverse terms
remain. It proves no global endpoint order identity can be inherited.
It does NOT establish its amplitude from the original seed or a 1%
failure in the required window. Initial append itself preserves colors.

## 4. What would settle the unchanged-law task

The EXACT missing quantities are the prepared-state vectors chi_+(t),
chi_-(t), their Gram entries, and the two-sector retention amplitudes,
uniformly over the frozen windows and stated actual record/receiver class.
For arbitrary inputs the worst output probability is min(p_+,p_-).
Conditional reference-entangled channel error additionally requires
quantitative control of their mismatch and Gram entries; exact equality
is required only for exact identity. Probability alone does not establish
the 1% task. No estimate for these growing-geometry quantities is claimed
here.

One sufficient route is still a uniform full-state comparison to the
endpoint reference, not the unsuccessful absolute majorant. If its error
is epsilon, the archived amplitude margins give

    output > (717/2000-epsilon)^2,
    retention > (659/1000-epsilon)^2,
    conditional trace distance <=epsilon/(717/2000-epsilon).   (4)

Output has the tighter PROBABILITY margin than retention. When a single
state-error bound is used for everything, the 1% INFORMATION condition
is more restrictive: epsilon<=717/202000 is sufficient for (4), with
both probability thresholds also satisfied. The parent target 1/1000
is safely below this ceiling. It remains UNACHIEVED, not a new pass gate.

At epsilon=1/1000 equation (4) yields output >0.12780625, retention
>0.432964 and trace distance <=2/715. These are conditional arithmetic
implications only. A direct bound of the complete instrument (2) could
succeed even if this stronger complete-state comparison is unsuitable.

## Verdict and limit of this cycle

There is no invariant fixed-record endpoint embedding of the supplied
two-fuel trajectory. The actual law couples it to an interior-source
environment. This closes the exact-invariance question, not the physical
operating question. The absolute Duhamel route cannot furnish the needed
full-window margin; the exact memory kernel and output Gram entries are
not controlled by the local identities. Finite primitive checks do not
enclose their long-time evolution, so no large-time conclusion is inferred.

The remaining problem is precisely (1)-(3) on the actual preparation and
frozen observation intervals. No source-only scalar theorem, unmeasured
receiver replacement, appended buffer or modified Hamiltonian resolves it.
No external review or physical validation occurred.

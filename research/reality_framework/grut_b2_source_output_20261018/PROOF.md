# B2-SO: a uniform original-source comparison

Author-derived; independent mathematical review pending. This is the
unchanged fixed-path matter law. It does not prove the full operating gate.

## 1. Physical task and one-source reference

Retain exactly CHARTER.md's H_L, double-fuel seed F_L, R_L, I_L and q=2
retention event. No conversion, return, charge state or collision term is
removed from H_L. Scalar notation suppresses ordered color, which factors
as I_4 for fixed endpoint fuels, including arbitrary entanglement with an
external reference.

Define the EXACT half-line one-source Hamiltonian h on fuel f and every
charged configuration |i,j>, i,j>=0. Its matrix is

    <i,j|h|i,j>=2+|i-j|+(1/4)1_(|i-j|=1)-A 1_(i=j=0),
    h_f,f=Omega, h_f,(0,0)=-J,
    plus hops=-j+, minus hops=-j-.

All reverse transitions and separations are included. One resource fixes
the capacity-one field. Set Pi=1_(-infinity,5/2)(h) and

    u(t)=exp(-ith)Pi f,    w=||(1-Pi)f||^2.             (1)

u is UNNORMALIZED. No physical low-band measurement or favorable prepared
packet is supplied; the difference from the original fuel seed is paid.

For even L, M=L/2-2, let C_M retain f and |i,j> with max(i,j)<=M. Define

    a_L(t)=||sum_(x=L/8)^(3L/8-1)|x,x><x,x| u(t)||^2,
    m_M(t)=||P_matter C_M u(t)||^2,
    b_M(t)=||1_(max(i,j) in {M,M+1}) u(t)||.            (2)

The receiver in a_L is INSIDE the cut for L>=24. Both a_L and m_M are
physical-coordinate probabilities of the exact filtered fuel evolution;
they include fuel/source interference and dressing, without renormalization.

## 2. Exact weighted invariant graph, not a finite-order dynamics

Let P span f and all co-pairs |n,n>, Q every separated charge state. Write

    h-2 = [ A_P  B* ; B  D+T ],
    D>=Delta=5/4, ||B||,||T||<=tau=5/96.

A_P has fuel diagonal -2A, root co-pair diagonal -A, conversion -J, and
zero other co-pair diagonals. Assign coordinate s_f=-1, s_(n,n)=n,
s_(i,j)=max(i,j) on Q. Primitive hopping and conversion change s by at
most one. For mu=3/10 use the symmetric weighted Schur norm

    ||K||_mu=max(sup_x sum_y |K_xy| exp[mu|s_x-s_y|/2],
                 sup_y sum_x |K_xy| exp[mu|s_x-s_y|/2]). (3)

This norm is submultiplicative, invariant under adjoint, and dominates
||exp(mu s/2)K exp(-mu s/2)||_2. It also applies to rectangular P/Q maps.
The associated matrix-kernel space is complete. Since

    exp(3/20)<=1+(3/20)/(1-3/40)=43/37<7/6,

we may use rational bounds

    ||B||_mu,||T||_mu <= b=35/576,
    ||A_P||_mu <= alpha=2A+7J/6.

Set x=1/19. The Riccati fixed-point map

    Phi(X)=-D^-1 B-D^-1 T X+D^-1 X(A_P+B*X)             (4)

maps the closed ||X||_mu<=x ball into itself because

    b/Delta+(b+alpha)x/Delta+b x^2/Delta < x.

Its Lipschitz bound is (b+alpha+2bx)/Delta<1/16. Thus it has a unique fixed
point in this ball. Let Z=[I;X], G=I+X*X and L_1=A_P+B*X. The fixed-point
equation gives

    h Z=Z(2+L_1),    ||L_1||_mu<=ell=alpha+bx<1/200.    (5)

D X=-B-TX+X L_1 is bounded, so the graph lies in the domain of h. The
closed graph is an invariant reducing subspace of the self-adjoint h.
On its orthogonal complement, vectors have form (-X*y,y), and

    <v,(h-2)v>/||v||^2
      >=(Delta-b-2bx-alpha x^2)/(1+x^2)>1/2.            (6)

The graph restriction is similar to bounded L_1 and has spectrum within
2+[-ell,ell]. Hence it is exactly the low spectral subspace Pi in (1),
not an assumed sector. Its orthogonal projector is Z G^-1 Z*. In
particular u(t)=Z exp[-it(2+L_1)]G^-1 f.

### A small correction bounds growth; it is not a phase approximation

Direct second-order algebra gives

    K_0=A_P-B*D^-1 B=-2A I-J(S+S*),                    (7)

where S shifts along f,0,1,... . There is only one physical boundary
counterterm, at event 0; no artificial end is added at M.

Use (4) once inside its T term. The term B*D^-1 T D^-1 B vanishes because
constituent hopping changes i+j parity and every co-pair has even parity.
The REMAINING exact correction E=L_1-K_0 satisfies

    ||E||_mu<=b x[(b/Delta)^2+(1+b/Delta)ell/Delta]<J/12. (8)

No term of E is discarded from the evolution in (5). Only its possible
contribution to weighted norm growth is bounded. Let W=exp(mu s/2).
The logarithmic norm of -i W K_0 W^-1 is <=2J sinh(3/20)<13J/42.
Including all of E, the exact weighted evolution obeys

    ||W exp(-it L_1)W^-1|| <= exp(c|t|), c=11J/28.     (9)

Weighted inverse control gives ||W G^-1 f||<=1/(1-x^2), and

    ||W u(t)||<=sqrt(1+x^2)/(1-x^2) exp(c|t|)
                 <(101/100)exp(c|t|).                (10)

Equations (4)-(10) retain the entire infinite charged space. They do not
exponentiate an energy-dependent Schur complement or accumulate a discarded
Hamiltonian error as |t| grows. Weighted norms can be justified first on
finitely weighted vectors and then by limits in the complete kernel space.

Therefore the actual filtered fuel preparation has the maximal-position
tail and boundary-shell estimate

    ||1_(max(i,j)>=M)u(t)||,
    b_M(t) <= (101/100) exp[-3M/20+c|t|].             (11)

Its cone coefficient c/(mu/2)=55J/21 is less than 3J, sufficient for the
pre-encounter cut. This is a statement about the PREPARED LOW PART; the
full fuel state also has high spectral components, paid once below.

## 3. High spectral weight is paid once

Q h Q>=h_* Q with h_*=2+Delta-tau. The P block is <=2+J. For real E
strictly between these edges, its Schur complement

    P h P-E-B*(Q h Q-E)^-1 B

is uniformly negative, so h has no spectrum there. This agrees with the
graph-defined low band. Every high spectral value is >=h_*.
The fuel seed has mean Omega and variance J^2. Its spectral measure gives

    w<=J^2/(h_*-Omega)^2=36/27613^2,
    sqrt(w)<=6/27613.                                (12)

No finite-relative-distance cutoff, bath reset or independent-emission
assumption is used. Filtered fuel and bulk co-pair preparation are distinct.

## 4. Joint comparison on the SAME physical Hilbert space

There is an isometry E_L from two cut one-source spaces into the full
two-resource space: keep the left coordinates 0..M, reflect the right
coordinates to L-1-M..L-1, and take the union of their constituent sets.
Map their fuel bits to the actual left/right fuels. Each factor contains
either one neutral charged pair or fuel. The intervening two events are
empty; each string is neutral outside its own interval. Capacity, Gauss,
hard-core conditions and ordered color are therefore respected exactly.

Set v(t)=C_M u(t), Phi_L(t)=E_L[v(t) tensor v(t)]. This is a mathematical,
unnormalized reference, NOT a second Hamiltonian or an implemented cut.

The compression E_L* H_L E_L equals h_M tensor I+I tensor h_M, where
h_M=C_M h C_M with its INHERITED diagonal. Cross-contact and field terms
vanish on this separated range. Do NOT use a newly terminated B1 path:
it would insert an extra -A co-pair energy at the cut.

The exact derivative defect has two sources: incoming amplitudes across
M+1 to M in each factor, and physical H_L hopping out of each cut. Each
boundary constituent shift has norm j+ or j-. Since ||v||<=1,

    ||(i partial_t-H_L)Phi_L(t)||
      <=2(j++j-)[||1_(max=M)u||+||1_(max=M+1)u||]
      <=4(j++j-) b_M(t).                            (13)

All non-boundary terms cancel as operators, not as classical histories.
Full evolution on the left side retains every subsequent collision,
reabsorption, two-source interference and coherent return.

Since C_M f=f, ||f-C_M Pi f||<=sqrt(w), and the product initial discrepancy
is <=2sqrt(w). Exact unitary Duhamel on the full H_L therefore proves

    ||exp(-itH_L)F_L-Phi_L(t)||
      <=delta_L(t):=2sqrt(w)+4(j++j-) integral_0^t b_M(s)ds. (14)

For 0<=t<=T_L=L/(6J), (11) and M=L/2-2 imply

    delta_L(t)<=12/27613+1500 exp(-L/105).             (15)

Indeed -3M/20+cT_L=3/10-L/105, exp(3/10)<10/7, and
4(j++j-)(101/100)(1/c)(10/7)<1500. The frozen sufficient threshold L>=2016
gives delta_L<1/1000. An exact Taylor lower bound on exp(96/5) verifies
this scalar comparison. Equation (15) holds on the ENTIRE [0,T_L], not
just one selected time. High-band norm is not charged per time step.

This is a uniform comparison theorem for the original source. It is not
an assertion that the reference's outgoing probability is large.

## 5. Linked bounds for the actual joint event

For L in 24N, the left receiver lies within 0..M. The reflected right
receiver is exactly the frozen right interval. The reference has

    ||R_L Phi_L(t)||=a_L(t),
    ||P_(q=2) Phi_L(t)||=m_M(t).

Thus (14), with NO multiplication of unproved physical successes, gives

    p_out^(L)(t) >= [a_L(t)-delta_L(t)]_+^2,
    p_q=2^(L)(t) >= [m_M(t)-delta_L(t)]_+^2.          (16)

The same physical R_L has the exact successful channel p_out rho on the
two ordered color qubits, including an external reference; failure is
erasure. No triplet substitution or fuel/global decoder is introduced.

The residual one-source quantities can be made still more explicit. Let

    g(t)=<f,exp(-ith)Pi f>.

Then, since all states except fuel contain matter,

    m_M(t)=1-w-|g(t)|^2-||(1-C_M)u(t)||^2.            (17)

For L>=2016 the last norm is <1/1000 throughout [0,T_L], by (11). A
SUFFICIENT pair of still-unproved original-source estimates is

    inf_(t in I_L) a_L(t)>=1/3,
    sup_(2/J<=t<=T_L) |g(t)|<=2/3.                    (18)

They would imply inf_I p_out>1/10 and inf_[2/J,T_L]p_q=2>1/4 after every
stated error is paid. Their constants are a sufficient proof interface,
not changed physical acceptance thresholds or observed probabilities.

## 6. A short launch corollary, not the growing-window result

B1's author-derived half-line fuel launch bound gives matter probability
>=22839/40000 for 2<=Jt<=4, under this same nominal one-source law. It is
a return-inclusive short-time statement, not its extrapolation to late
times. Subtract the one-source high norm <=6/27613<1/4000 and the tail in
(11). For L>=2016 and Jt<=4 the latter is <10^-6. It follows that

    m_M(t)>57/100,
    p_q=2^(L)(t)>(57/100-1/1000)^2=0.323761,
                                      2<=Jt<=4.     (19)

The rational root checks are (3778/5000)^2<22839/40000 and
(3778/5000-1/4000-10^-6)^2>57/100. This corollary depends on B1's prior
proof as well as the present author-derived comparison. It does not put
matter in the distant L-scaled receivers or cover retention after theta=4.

## 7. Unresolved estimate and stop

This cycle proves the proposed joint comparison, subject to independent
review, and pays its growing-time error. It does NOT prove either inequality
in (18). The exact fuel Green function is

    G_f(z)=1/[z-Omega-J^2 <0,0|(z-H_ch)^-1|0,0>].      (20)

Its source spectral measure and outgoing generalized eigenfunctions must
bound the receiver probability a_L and fuel return g. Available high-band
weight and bulk group velocities alone do not provide these lower/upper
bounds. Applying only positivity and norm conservation to (2) permits
a_L=0 and yields a zero operating lower bound in (16).

A rigorous boundary spectral-density/scattering estimate, including
localized components, slow modes, physical co-pair detection and finite-time
remainders uniformly over I_L, is the remaining input. The source-to-output
lemma therefore stays unresolved, and the full operating gate is NOT passed.
No empirical novelty, extensive medium, moving-geometry theorem, parameter
robustness or successful repaired B2 is claimed. Stop this cycle here.

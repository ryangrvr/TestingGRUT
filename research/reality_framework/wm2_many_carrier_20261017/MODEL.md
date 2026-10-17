# WM2-D1: explicit finite-density extension and its capacity price

**New supplied candidate law; author-derived analysis, independent review
pending.** D1 is not unchanged UC1/IG1. Relational locality remains a primitive.

## 1. Why the literal registers cannot support the requested limit

UC1 fixes N_A,N_B<=1 and f_A,f_B in {0,1}, hence P=N_A+N_B<=2. Even
enlarging only those counters is insufficient. On a rooted occupied path,
Gauss and flux support give, in root-to-head orientation,

    E_(x,x+1) = -sum_(y>x) (n_Ay+2n_By).                       (M1)

All charges are positive. A link capacity |E|<=C therefore forces the total
charge away from the root to be <=C, by the first edge. For hard-core slots
with at most one of each species at an event,

    P_away<=C, P_total<=C+2.                                   (M2)

At UC1's C=3 this is at most five particles after relaxing its original
two-particle cap. No extended charged density follows at increasing N.
M2 is specific to positive charges, a simple rooted path and compensation
at its root; it is not a no-go for distributed compensators or neutral matter.

## 2. D1 registers and exact physical space

At every event supply a hard-core A slot with states {vac,0,1} and a hard-core
B slot with states {vac,0,1,2}. These are commuting local tensor factors,
dimension 12 per event, not a derivation of fermionic statistics. They permit
at most one of each species per event, two carriers total there.

Let N_s=sum_i n_si and enlarge f_s to the integer range 0,...,N. Enlarge
every potential unordered-pair electric register to E=-3N,...,3N. Keep
U_(i->j) a partial shift lowering the oriented E by one, with its exact
adjoint; [E,U]=-U also at the boundary. Retain the commuting constraints

    G_i=sum_(j!=i) E_ij -n_Ai-2n_Bi+R_i(f_A+2f_B)=0,
    C_ij=(1-B_ij) E_ij^2=0, N_s-f_s=0.                         (M3)

B_ij is the current occupied consecutive bond, not supplied label adjacency.
The physical state space is the joint kernel. Every positive carrier must
lie on the current word. Its fields are uniquely M1, with magnitude at most
3(N-1). Thus capacity 3N covers the **entire** declared hard-core physical
space, with all number sectors and internal states, not a chosen low-energy
packet. No physical electric configurations have been truncated numerically.

The enlargement costs O(N) electric levels per potential link (O(log N)
qubits per link) and O(N) compensator levels. There are still O(N^2) potential
pair registers. Shift operator norms remain one; their dimension growth does
not increase the explicit physical hopping matrix element. The root charge
f_A+2f_B can be extensive. None of these resources is dynamically selected.

## 3. Hamiltonian — all reverse processes retained

Retain the WM2 storage/receiver law, L_N=64N^2 modes, receiver charge
n_R=|q|-1, delta_i=bar_delta+epsilon_i with |bar_delta|<=1/2,
|epsilon_i|<=1/4, and common 0.9<=g<=1.1. Keep zero electric energy.
With epsilon=1/4 retain every zero-flux-gated assisted append/pop (g/N)
and bare append/pop (epsilon/sqrt(N)), using the head/free or head/parent
gate n_A+n_B. This is a hard-core occupation gate, bounded by two even
when P is extensive. Keep B's on-event level-2 clock, with coefficient one.

Let b_sia annihilate internal level a in species s's slot, z_A=1,z_B=2.
The explicit matter witness is

    H_hop=-sum_(occupied {i,j},s,a) j_s
            (b_sja^* U_(i->j)^(z_s) b_sia + adjoint),
    j_A=1/64, j_B=1/96,
    H_bond=(1/4) sum_(i!=j) B_ij Z_Ai X_Bj.                    (M4)

Z_A and X_B are zero in their slot vacuum; X_B acts only on levels 0,1.
The third B level is present and mobile. Co-located species are permitted.
Each term has an explicit ambient-register implementation preserving M3.
A pop through a charged head is blocked by zero flux; no hidden P h P gate
or constraint-violating virtual state implements the Hamiltonian.

For a general finite-range class replace M4 by explicitly implemented,
number-conserving local matter generators h_I on occupied intervals I of
diameter <=R, with the IG1 restrictions on coefficient-reading support and
no new direct receiver/storage channel. On each fixed word their exact tree
isometry must give ordinary local matter operators, and h_I annihilates its
interval vacuum. Replace IG1's global w cap by

    sup_word sup_x sum_(I contains x) ||h_I|| <= J_0.           (M5)

Norms in M5 are on the **full local hard-core space**, not just a favorable
number sector. On-site position-preserving internal energies can be treated
exactly; the primary class retains the unit B clock. The primary theorem uses
R=1,J_0=2. These are supplied bounds on interactions, not consequences of G.

For M4 a bond has norm <=j_A+j_B+2kappa=101/192, kappa=1/4. Therefore
J_0<=101/96<2 on every word and every reflecting prefix. The rule admits
interacting internal states and both charged species at arbitrary hard-core
occupation. It does not imply equal speeds or nonzero response in jammed states.

## 4. Extensive norm is real — and not a propagation velocity

The primary reference instance fixes g=1 and delta_i=0 for every event;
all other coefficients are the fixed values above and all N-dependences are
the declared formulas. The independent-disorder interval is a robustness
class, not permission to fit a different assignment after each outcome.
Finite diagnostics deliberately use g=1 and delta_i=.23 plus equally spaced
offsets in [-.25,.25], an assignment inside that frozen class. The free-channel
holdouts are uniform in the class, so neither g nor delta is fitted to them.

On fixed N_A,N_B, P=N_A+N_B, the changing witness matter norm obeys

    ||V_mat|| <=2j_A N_A+2j_B N_B+2kappa min(N_A,N_B)=O(P).      (M6)

Hopping has at most two outgoing moves per carrier; internal color adds no
extra outgoing channel for an occupied slot. The internal interaction has
at most two occupied neighbors per carrier. Hermitian row sums give M6.

For N_A=N_B=n, place alternating A,B on 2n consecutive sites, A internally
Z=+1 and B internally X=+1. Hopping has zero expectation in the fixed-position
state and the bond expectation is kappa(2n-1). Hence ||V_mat||>=kappa(P-1).
The global norm is Theta(P), including on sectors P proportional to N.
Using it in IG1's total-norm Dyson estimate would shrink a length-N window
to order N/P. D1 instead uses M5 and connected propagation estimates.

The geometry-changing norm is still <=g+4epsilon, independent of P:
the directed assisted map has squared column norm <=(g/N)^2(N-|q|)|q|
and unique parents; the directed bare map is bounded by 2epsilon because
the head occupation is at most two. Adding adjoints gives g+4epsilon.
Unbounded occupation statistics would require a different bound.

## 4a. Self-contained operator and preparation dictionary

The geometry basis contains every simple rooted ordered word w, with no
repeated label, 1<=m=|w|<=N. Let a_i(w)=1_(i in w). Storage q_i is a bit,
q_root=1, and the full allowed bosonic receiver occupations obey n_R=|q|-1;
in particular none is discarded as a 'storage error'. With h_L the open
L-site matrix having off-diagonal -1, the relative dimensionless diagonal is

    D_0=16(|q|-m)+64 sum_i(q_i-a_i(w))^2
          +sum_i delta_i(q_i-R_i)+dGamma(h_L),
    D=D_0+sum_i n_Bi |2><2|.                                 (M7)

At a current head h and free event v, the assisted forward primitive appends
v, raises q_v from 0 to 1, creates a contact boson with amplitude
-(g/N)sqrt(n_0+1), and carries Pi_(E_hv=0). The bare forward primitive appends
v, leaves q and receiver unchanged, and has coefficient
-(epsilon/sqrt(N))(n_h+n_v) Pi_(E_hv=0), n_i=n_Ai+n_Bi. Both use the original
head/free/bond register operations, acting on h,v and receiver contact as
declared, not a supplied label ordering. Their exact adjoints are the pops.
On the physical kernel n_v=0 at a free event. These primitives plus M4 and
M7 specify H_D1 completely; general M5 gives a declared wider class, not an
unfixed term in the concrete witness.

This is a relative energy convention. Replacing h_L by 4+h_L and adding
-4(|q|-1) to the storage diagonal changes nothing on the exact charge sector
n_R=|q|-1 and gives receiver frequencies in [2,6]. No implicit energy source
or irreversible sink is inserted by that rewrite. Subsequent closed evolution
conserves the specified total H_D1 energy.

The neutral preparation is the original uniform-root bondless one-head state,
q_i=R_i, receiver vacuum, all matter slots vacuum, f_s=0 and E=0. Root components
are orthogonal; every stated uniform formation bound is root-uniform and may
therefore be averaged with those original weights. An optional inert CP-12
benchmark factor is a spectator, not an energy supply to D1. The populated
preparations and pulse/readout observables are specified in PROOF.md; the
fixed empirical-channel holdouts are in EFFECTIVE_PREDICTIONS.md.

## 5. Formation and preparation are different statements

In joint matter vacuum, f_s=0 and all fields vanish. Every new matter term
and bare occupation gate is zero. The old WM2 formation law is exactly
recovered, with no change to S9's author-derived, independent-review-pending
status or large receiver/time resources.

But [H,N_s]=0. The original neutral seed has **zero probability of creating
P>0** under this law. A hard-core seed with P extensive already requires
at least P/2 occupied events by M3. D1 does not autonomously generate a
populated macroscopic phase from the original one-event neutral seed.

The positive regime in PROOF.md starts after recording and after explicit
matter/countercharge/flux preparation. It retains the actual conditional
receiver state. Neither the success probability nor the always-on implementation
of that extensive preparation is derived. It is not an original-trial formation
and matter-creation probability. The extension answers a conditional finite-
density use question, not that separate formation problem.

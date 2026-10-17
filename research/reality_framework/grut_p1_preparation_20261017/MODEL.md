# P1 microscopic registers and controlled Hamiltonian

**Separately supplied model; author-derived analysis, independent review
pending.** The primary architecture is a post-recording, externally powered
protocol. Locality of its interaction generators is supplied. It is not an
autonomous or parameter-free theory.

## Registers and physical space

Retain every simple rooted word w of length m=1,...,N, the original head,
free/bond flags, storage bits q_i with q_root=1, and all receiver occupations
with n_R=Q-1, Q=sum q_i. Put a_i=1_(i in w), R_i=1_(i is root). The physical
receiver has L=64N^2 modes, open-chain one-particle hopping h_L with adjacent
entries -1. All storage-inconsistent words are retained.

Replace D1's positive species/countercharge registers by two commuting
hard-core slots at each event, c_+,c_-, each with states {vac,occupied}.
Co-location is allowed; n_i=n_+i+n_-i<=2 and charge rho_i=n_+i-n_-i.
No compensating root register remains. This does not derive fermion statistics
or the observed particle spectrum.

Every potential unordered pair has a three-level electric register. In
orientation i->j, E_ji=-E_ij, E_ij in {-1,0,1}; U_(i->j) lowers E_ij by one,
annihilating its lower boundary. Its adjoint raises it. In particular,
[E_ij,U_(i->j)]=-U_(i->j), including both boundaries. A label ordering may
store this orientation in code but does not occur as a physical affinity.

Every potential unordered pair also has a fuel qubit f_e, n_fe=f_e^*f_e,
energy Omega_e=3+Delta_e with |Delta_e|<=1/2. Inactive fuel qubits start in
vacuum and remain there. They need not lie on an occupied link.

Use the exact joint kernel of

    G_i=sum_(j!=i)E_ij-rho_i,
    C_ij=(1-B_ij)E_ij^2,
    M_i=(1-a_i)n_i,
    S=Q-1-n_R.                                                (1)

M_i is a **new supplied matter-support constraint**. With opposite charges,
Gauss and flux support alone would allow a neutral co-located pair on a free
event. The empty-head gate below is essential, not an inherited D1 fact.
Summing G_i gives N_+=N_- on every physical state. On a path, the unique field
is E_(x,x+1)=sum_(y<=x)rho_y=-sum_(y>x)rho_y. Keep every matter configuration
satisfying this relation and |E|<=1, not just adjacent dipoles. The fixed
capacity is part of the new law, not a numerical truncation of D1.

## Diagonal and geometry terms

At all stages retain

    D_geo=16(Q-m)+64 sum_i(q_i-a_i)^2
                   +sum_i delta_i(q_i-R_i)+dGamma(h_L),
    H_E=sum_(unordered e)E_e^2,
    H_mass=N_++N_-, H_f=sum_e Omega_e n_fe.                     (2)

Here delta_i=bar_delta+epsilon_i, |bar_delta|<=1/2,
|epsilon_i|<=1/4. The electric coefficient is fixed and positive: one.
The positive-frequency receiver convention replaces h_L by 4+h_L and adds
-4(Q-1) to storage, exactly the same Hamiltonian since S=0. Receiver modes
then have energies in [2,6]. This is not a new source of fuel energy.

Assisted append at head h and free event v raises q_v, creates one receiver
contact boson, and appends v, with coefficient -g/N, g in [0.9,1.1]. Add
Pi_(E_hv=0) and Pi_(n_v=0) to the forward primitive and its exact adjoint.
Thus an assisted pop requires its removed head to be empty as well as its
removed link to have zero field. Bare append/pop has coefficient
-(epsilon/sqrt(N))(n_h+n_v), epsilon=1/4, the same electric/empty-child gates,
and leaves q/receiver unchanged. The forward empty-child gate becomes the
empty-head gate in the adjoint. All reverse transitions are retained.

The norm of the entire geometry-changing operator is <=g+4epsilon<=2.1.
For assisted creation the squared directed column norm is at most
(g/N)^2(N-Q)Q<=g^2/4, and distinct parents have disjoint children. Empty-child
gates only reduce it. For bare creation n_h<=2 and there are <=N children;
directed norm <=2epsilon and adjoint doubling gives 4epsilon. No extensive
matter norm enters this estimate.

In matter and fuel vacuum, all new diagonals, matter terms and bare gates
vanish, and the empty-child gate equals one. The formation embedding is
exactly the archived WM2 vacuum law. It inherits S9's author-derived status,
not independent validation.

The complete original auxiliary seed is
N^(-1/2) sum_r |w=(r),q_i=1_(i=r),receiver vacuum,matter vacuum,E=0,fuel vacuum>.
Root components are orthogonal and retain these original uniform weights;
all conditional bounds are root-uniform before averaging. An optional old
CP-12 benchmark factor is an inert tensor factor, not a source of energy
or a claimed common matter sector. The finite code uses the root-0 component
only, with relabeling handled by the analytic uniformity, not numerical sampling.

## One fixed control schedule

At T_f=4N^2 measure the word, retaining the actual conditional receiver
state. Accept m>=N/2; rejected trials contribute zero. Define K=N/32 and
the matching M={(w_0,w_1),...,(w_(2K-2),w_(2K-1))}. The controller reads
recorded rank, not a predetermined event-label ordering. This is supplied
feed-forward and requires a record/addressing system; it is not autonomous
selection of an interaction law.

For local stage time s in [0,pi/2], use

    H_L(s)=D_geo+H_geo+H_E+H_mass+H_f
              +sum_(e in M)eta_e(e^(-3is)f_e^*+e^(3is)f_e),    (3)

where eta_e in [0.9,1.1]. This fixed-frequency drive is the external work
source. Initially every fuel is in vacuum. The phase convention is reset
by the declared controller clock at this stage's beginning; it is not
adjusted to Delta_e. Geometry remains on.

For the next pi/2, turn off that drive and use

    H_C=D_geo+H_geo+H_E+H_mass+H_f+V_C,
    V_C=sum_(e=(i,j) in M)lambda_e B_(i->j)
          (f_e U_(i->j)^* c_+i^* c_-j^* + adjoint),            (4)

lambda_e in [0.9,1.1]. B_(i->j) is the current directed occupied bond;
the selected orientation is its recorded root-to-head orientation. It
preserves the two charges' opposite signs and raises E_ij by one. It
also permits recombination and exact fuel return. Both directions carry
the same pulse envelope. No independent emission assumption is used.

At the end of this pulse switch V_C off in both directions, and switch on

    V_op=-sum_(occupied {i,j},sigma)j_sigma
              (c_sigma,j^* U_(i->j)^(z_sigma)c_sigma,i+adjoint)
           +kappa sum_(occupied {i,j})
                 (n_+i n_-j+n_-i n_+j),                     (5)

where z_+=1,z_-=-1, U^(-1) means U^*, j_+=1/64, j_-=1/96,
kappa=1/4. The first sum uses one orientation per bond plus its adjoint.
The operating Hamiltonian is the base of (3) plus (5). N_+,N_- are now
conserved. Persistent matter number relies on this explicit switch;
P1 does not establish persistence with conversion always on. Recombination
was included throughout the conversion interval, not replaced by a rate.

All stages and switches are part of the declared protocol. No geometry
term is switched off or clamped. No bath, receiver reset, matter insertion,
or measurement selecting a favorable empty buffer is used.

## Constraint checks as operators

For a positive hop i->j, rho_i decreases and rho_j increases, while U lowers
E_ij; the change of sum E equals the change of rho at both endpoints.
A negative hop uses U^*, reversing both signs. For conversion, creation
changes (rho_i,rho_j) by (+1,-1) and U^* changes the two divergences by
(+1,-1). Therefore [G_i,H]=0 for every stage, including the drive.

Matter hops and conversion contain their occupied-bond projector, hence
commute with C and M on the ambient valid-word register space. Geometry
changes neither charge nor field; its zero-field gate makes C commute
through the changed bond and its empty-child/head gate makes M commute
through the changed membership bit. Assisted geometry changes Q and n_R
together; other terms change neither. Thus the kernel of (1) is invariant
under every stage and switch. These are explicit ambient-register operations,
not arbitrary P h P implementations. The finite shift is not a cyclic link.

During conversion F+N_+=F+N_- is conserved, F=sum n_fe. Loading keeps
N_+=N_-=0 and F<=K. Consequently all subsequent physically reachable
states have N_+=N_-<=K and F+N_+<=K. All number sectors, fuel coherences,
geometry excursions, storage errors and receiver returns within those exact
invariants remain in the analysis.

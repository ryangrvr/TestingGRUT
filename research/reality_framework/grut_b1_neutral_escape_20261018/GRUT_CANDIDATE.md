# GRUT-B1: one complete conditional candidate

Author-derived, independent proof review pending. Quantum mechanics, external
time, the path assembly grammar, and universal occupied-interval support of
interaction generators are supplied. No derived universal locality is claimed.

## Registers and physical space

Retain P1's rooted simple word w, attached bits a_i, storage bits q_i,
receiver with L=64N^2 physical modes, and root flag R_i. All storage-error
states and reverse append/pop transitions remain in the law. On every
potential event pair, E_ij=-E_ji has capacity -1,0,1; B_ij is occupied-bond
adjacency. On every event put a hard-core plus with states 0,up,down, a
hard-core minus with states 0,1, and a neutral fuel with states 0,up,down.
f_i,sigma=|0><sigma|. Each species has at most one particle per event;
opposite charges may co-locate. Distinct events/species commute.

Constraints, with n_+ summed over its two colors:

    G_i=sum_j E_ij-n_+,i+n_-,i=0;
    (1-B_ij)E_ij^2=0;
    (1-a_i)(n_+,i+n_-,i)=0;
    (1-a_i)n_f,i=0.

The last fuel-support condition is new; all append/pop empty-head gates
include fuel as well as matter. No evolution may orphan a fuel on a free
event. The joint kernel is invariant. Let U_(i->j) lower E_ij by one,
annihilating a capacity violation. These partial shifts obey [E,U]=-U.

## Fixed Hamiltonian, all times

H_B1=H_geo + H_E + H_mass + H_density + H_hop + H_fuel
     +H_degree + H_conversion.

* H_geo is the S9/P1 storage/receiver law with P1's empty/zero-field gates,
  strengthened by the fuel empty-head gate. In relative-energy convention,
  its diagonal is 16(Q-m)+64 sum_i(q_i-a_i)^2
  +sum_i delta_i(q_i-R_i)+dGamma(h_L), Q=sum_i q_i,
  n_receiver=Q-1. h_L has nearest-mode hopping -1. Assisted append creates
  one receiver contact boson and raises q_child, coefficient -g/N;
  bare matter-driven append has coefficient
  -(n_head+n_child)/(4 sqrt(N)); exact adjoints are retained. The old
  storage-error sectors are not projected out. |delta_i|<=3/4, g in [.9,1.1]
  with the inherited common-plus-independent disorder prescription.
* H_E=sum_occupied E_e^2; H_mass=sum_i(n_+,i+n_-,i).
* H_density=(1/4)sum_occupied(i,j)(n_+,i n_-,j+n_-,i n_+,j).
* H_hop=-sum_occupied(i,j),sigma j_+
  (c^dagger_+,sigma,j U_(i->j)c_+,sigma,i+h.c.)
  -sum_occupied(i,j) j_-(c^dagger_-,j U^dagger_(i->j)c_-,i+h.c.).
  j_+=1/64, j_-=1/96.
* H_fuel=Omega sum_i n_f,i.
* H_degree=-A sum_i(2-deg_w(i)) n_+,i n_-,i,
  A=(j_+^2+j_-^2)/(5/4)=13/46080.
* H_conversion=-lambda sum_i,sigma
  (f_i,sigma c^dagger_+,sigma,i c^dagger_-,i+h.c.).

Only a supplied fuel-site program has conversion enabled. Primary program:
the root alone. Secondary two-cell diagnostic: both ends of the short path.
This program is static, not a clock. The on-event neutral conversion is a
new primitive operator, not P1's adjacent-charge conversion. It commutes
with Gauss without an electric shift. Both directions remain on forever.
All terms preserve the constraints. The degree counterterm, color-blind
kinetic terms, fuel energy and enabled site are explicit model assumptions.

J=2j_+j_-/(5/4)=1/3840; Omega0=2-2A.
Primary class: |Omega-Omega0|<=J/1000, |lambda-J|<=J/1000.
This is supplied renormalized resonance over a narrow size-independent
interval; no principle selecting it is derived.

## Preparation and resources

Primary conditional experiment: accept a geometry record with m>=N/2,
N>=2^28, retaining its actual conditional receiver state and q_i=a_i.
Set matter and E to vacuum, prepare one excited root fuel with arbitrary
unknown internal density matrix rho. Then evolve this one fixed Hamiltonian.
The record, fuel excitation and internal input are supplied interventions;
there is no later switch or reset. No desired matter state is inserted.

Before fuel loading, matter/fuel vacuum is an exact invariant block equal to
the old vacuum formation law. Thus S9 can provide a preparation corollary
*if its proof survives review*. It does not prove fuel loading or recording
autonomous. Geometry remains fully quantum during use; only a mathematical
comparison retains a prefix of ell=N/4 events.

Exact conserved resources:

    N_+=N_-=1-F,  F=sum_i n_f,i in {0,1};
    total net charge=0.

The input qubit resides either in the fuel or the one plus. Color-blind
evolution factorizes exactly as H_spatial tensor I_2 in this sector. This
statement is not inherited at finite density with indistinguishable carriers.

Initial added mean energy is Omega, variance lambda^2; mass/field/matter
interaction expectations initially vanish. It is not an energy eigenstate.
Fuel energy is the work resource, and its energy spread is stated. Total
energy after preparation is conserved. Electric energy is positive and at
most m-1 on a recorded path; in the one-resource, low-energy experiment
its expectation is O(1), not a root-compensated extensive field. The small
negative degree term is bounded below by -2A per pair and cannot cause an
extensive negative-energy instability. The inherited receiver energy
convention is (4+h_L) with compensating -4(Q-1) storage energy, not a second
fuel supply. Full potential-register count is still O(N^2), receiver L=64N^2.

## Measurements and claimed domain

At theta=Jt=4, report successful arrival only if the retained recorded prefix
survives and a co-located pair is at a recorded position x>=4 from the root.
Decode the plus color there. Trace other degrees and flag every failure as
erasure. No favorable receiver reset or arrival renormalization. In this
one-resource sector the channel is p_arr rho plus an orthogonal erasure of
weight 1-p_arr, with unit conditional fidelity for arbitrary qubits.

The claimed operating interval is theta in [2,4]: t in [7680,15360].
This is a conditional one-composite quantum medium. No extensive matter,
collision lifetime, continuum spacetime, gravity, or distinctive empirical
prediction is claimed. Finite-path recurrences are permitted outside it.

# A1 — complete fixed law and initial resource

A1 is a minimal successor to controlled P1, not a new interpretation of P1.
It starts from a supplied record and fuel, not the original empty seed. All
claims here are author-derived, independent mathematical review pending.

## Registers, constraints and terms

Keep P1's rooted-word states w, a_i(w), B_ij(w), root bits R_i, storage q_i,
Q=sum q_i, receiver n_R=Q-1 with L=64N^2 modes, and every storage-error state.
Two commuting hard-core slots c_+,c_- have charges +1,-1 and unit mass energy.
Every potential pair has E=-1,0,1 with partial shift U_(i->j) lowering E_ij,
E_ji=-E_ij. Every pair also has a fuel qubit f_e. Keep the physical kernel

    G_i=sum_j E_ij-n_+i+n_-i=0,
    C_ij=(1-B_ij)E_ij^2=0, M_i=(1-a_i)(n_+i+n_-i)=0,
    Q-1-n_R=0.                                                (1)

Add static program bits p_(i->j), which commute with H and are supplied in
the recorded matching described below. They are not another dynamic selector
or a clock. The constant law is, for all t>=0,

    H_A1=D_geo+H_geo+sum_e E_e^2+(N_++N_-)+sum_e Omega_e n_fe
          +V_hop+V_den+V_convert,                             (2)
    D_geo=16(Q-m)+64 sum_i(q_i-a_i)^2
                      +sum_i delta_i(q_i-R_i)+dGamma(h_L),
    V_hop=-sum_(occupied {i,j},sigma)j_sigma
                 (c_sigma,j^*U_(i->j)^(z_sigma)c_sigma,i+h.c.),
    V_den=kappa sum_(occupied {i,j})(n_+i n_-j+n_-i n_+j),
    V_convert=sum_(i->j)p_(i->j)lambda_(ij)B_(i->j)
                 (f_(ij)U_(i->j)^*c_+i^*c_-j^*+h.c.).

z_+=1,z_-=-1, U^(-1)=U^*, j_+=1/64,j_-=1/96,kappa=1/4.
Each occupied bond is counted once in the hopping sum, plus its adjoint.
There is no loading drive or time envelope. In particular recombination
and fuel return remain in (2) during matter transport.

H_geo is P1's assisted append/pop of coefficient -g/N, g in [0.9,1.1],
raising/lowering q and a receiver contact quantum together. Bare append/pop
has coefficient -(1/(4sqrt(N)))(n_h+n_v), leaves storage/receiver unchanged.
Both carry zero-field and empty-child gates in the forward append and their
exact adjoints; a pop requires its removed head to be empty. Their total
changing norm is <=g+1. No physical geometry state is discarded. Use
delta_i=bar_delta+epsilon_i, |bar_delta|<=1/2, |epsilon_i|<=1/4.

Every term is the explicit P1 ambient-register implementation; the constraint
commutators stay zero because the previously separate terms now are summed.
No new projection implements the dynamics. Partial shifts, all reverse
processes, neutral co-location, and all physical field sectors remain.

## Supplied initial condition

Take N a multiple of 64, K=N/64, a recorded length m>=N/2, q_i=a_i, and
its actual conditional receiver state, unreset. At recorded positions

    (4q,4q+1,4q+2,4q+3), q=0,...,K-1,

put the selected program edge between 4q+1 and 4q+2 and its fuel in |1>.
The two outer events are guards, initially empty. All other program/fuel
registers are zero. All matter and fields start vacuum. The program and
selected excited fuels are a supplied state, not prepared autonomously in
this campaign. A degenerate classical record/program needs O(N log N) bits;
addressed fuel count is O(N). No clock register is used.

Reference parameters are lambda_e=1 and Omega_e=13/4. The permitted class
also contains independent lambda_e in [0.9,1.1], Omega_e in [5/2,7/2].
Fuel energy 13/4 matches the adjacent-pair energy 2+1+kappa; the additional
kappa is essential when V_den remains on. This is a supplied resonant point
inside a nonshrinking class, not a selected natural constant. A counterexample
at this point suffices to reject a retention guarantee uniform on that class.

From this state the EXACT invariants are

    N_+=N_-, F+N_+=K, F=sum_(selected e)n_fe,
    P=N_++N_-=2(K-F), 0<=P<=2K.                              (3)

All unselected fuel registers stay vacuum because their program bits are zero.
The zero-fuel, zero-matter seed remains matter vacuum exactly. A1 therefore
does not claim to create matter without a supplied fuel resource.

## Energy and authority prices

The initial added mean energy is sum Omega_e=13K/4 at the reference point.
The conversion expectation is initially zero; its initial variance is
sum lambda_e^2=K. Orthogonal one-cell conversion images establish this
identity. The resource state is not an eigenstate of the already-on H:
its energy spread and preparation are explicit inputs. Geometry/receiver
energy is the inherited conditional account, not a second fuel supply.

After t=0 total H_A1 energy is conserved and no controller work, reset,
switch or measurement feedback occurs. Electric energy is bounded by N-1,
fuel by (7/2)K and mass by 2K on every reachable state; changing and density
norms are O(N), with fixed local coefficients. Positive receiver frequencies
4+h_L and storage offset -4(Q-1) give the same original relative Hamiltonian.
There is no root countercharge, new bath or implicit dissipative rate.

Original formation/recording/fuel preparation costs remain supplied. If one
wishes to connect to S9, those interventions and its independent-review
dependency must be restored explicitly. The retention counterexample below
does not require S9 and cannot be repaired by increasing its acceptance rate.

# Local formation and an extensive return obstruction

Author-derived; independent proof review pending. The obstruction uses the
permitted point lambda=1, Omega=13/4, retaining BOTH nonzero charged hoppings,
all nearest-bond density interactions, positive field energy and recombination.
It is not a no-go for every autonomous construction or every detuning.

## 1. One exact conversion cell

Ignoring only the hops for the moment, the initial fuel |F> and oriented
adjacent pair |D_+> have the same diagonal energy Omega=2+1+kappa=13/4.
Their Hermitian block is Omega I+lambda sigma_x. Pair probability is
sin^2(lambda t); it is one at pi/2 and zero at pi, lambda=1.

The actual two-event one-resource physical sector has FIVE states:
F, D_+=(+left,-right), D_-=(+right,-left), M_left and M_right
(co-located neutral pairs). Diagonal is (Omega,Omega,Omega,2,2).
Conversion connects F,D_+ with lambda. Hopping connects

    D_+ <-> M_left: -j_-,  D_+ <-> M_right: -j_+,
    D_- <-> M_left: -j_+,  D_- <-> M_right: -j_-.

This is an exact finite Hamiltonian, not a one-way rate or discarded channel.
Writing a=j_++j_-, b=j_+-j_-, its exact fuel Green function is

    G_F(z)=[z-Omega-lambda^2 g_D(z)]^(-1),
    g_D(z)=(z-2)/2 * {1/[(z-Omega)(z-2)-a^2]
                       +1/[(z-Omega)(z-2)-b^2]}.

Its poles/residues give the exact fuel amplitude sum_l r_l exp(-i E_l t).
No continuum or sink occurs in this cell. The bounded matrix controls check
this formula and its nonzero-hopping evolution directly.

## 2. A local isolation estimate without an extensive error

First work on the fixed retained prefix, keeping the links. For any selected
four-event cell, define H_cut by removing only the five hopping bonds that
touch its four events (fewer at a boundary). Keep every density term,
conversion term and diagonal, and let the rest of the world follow its full
interacting evolution. This is a mathematical reference only.

The cell's guard events then stay empty. Cross-boundary density products
annihilate them. Each separated outside component has zero conserved charge,
so the cut fields remain zero. Therefore the cell evolves as the exact
two-state conversion block, tensor its empty guards, even while the outside
contains arbitrary hopping, conversion, matter interactions and returns.
This uses its declared initially empty guards, not a later favorable outcome.

Let V_cut=H-H_cut. Its five consecutive hopping bonds lie on at most six
consecutive events. On a fixed tree, solving Gauss identifies its kinetic
matrix with the compression of ordinary two-species hard-core hopping to
|E|<=1 configurations. This is ONLY a norm estimate; positive E^2 is not
treated as a local matter-only operator, and no projected law is substituted.
For one species on six sites, Jordan-Wigner/free-path diagonalization gives
the many-body kinetic norm as the sum of its three positive one-particle
energies. Their squared sum is 5 j^2 (half trace of the one-particle square).
Cauchy-Schwarz bounds their sum by sqrt(15)j<4j. Compression cannot enlarge
that norm. Consequently, uniformly over the rest of the system,

    ||V_cut||<=4(j_++j_-)=5/48=:eta.                          (4)

Boundary cells with fewer bonds obey the same bound. The proof includes
the finite partial shifts and all allowed particle numbers. It does not
assume independent physical cells or use ||H||=O(N) as a small perturbation.

Exact unitary Duhamel gives ||(U-U_cut)psi||<=eta |t|.
At t_f=pi/2 the reference fuel is empty, so

    <n_fe(t_f)> <=(eta*pi/2)^2.

At t_r=pi the reference fuel is full, so

    <1-n_fe(t_r)> <=(eta*pi)^2.                              (5)

These follow from projection AMPLITUDE bounds before squaring. An expectation
norm bound alone would be weaker. Sum (5), use the exact resource identity
P=2(K-F), and apply Markov to F at t_f and P at t_r:

    Pr(P(t_f)>=K) >=1-eta^2*pi^2/2 >0.946,
    Pr(P(t_r)>=K) <=2 eta^2*pi^2 <0.215 <1/4.                 (6)

These size-uniform statements concern an extensive threshold. No physical
factorization, independent-emission process or thermodynamic limit was assumed.
They show formation followed by substantial loss under the SAME fixed law.

## 3. Complete geometry and receiver comparison

The word is not frozen physically. Retain ell=N/4, b=N/16-1 and a=m-ell>=N/4.
On the common physical Hilbert space use the same reflected-tail reference
method as P1, now including all conversion terms inside the prefix. Keep
the actual receiver and all storage/field states in the physical law.

The conserved local resource is different from P1's post-switch particle
number. Assign each selected fuel to its parent and set

    q_x=n_+x+n_-x+2n_f,x <=4,
    W_mu=sum_x exp(mu(x-b))q_x, C_mu=4/(1-exp(-mu)).

Every local conversion/hopping conserves sum_(x in its bond)q_x. Diagonal
terms commute with W_mu. Initially W_mu<=C_mu, even with K extensive.
With links retained, the ordinary bond norm is <=293/192; the additional
one conversion term at a participating site has norm <=1.1. Thus the local
interaction bound J_0=6, range one, is sufficient on the full local space.
The weighted-form/Gronwall argument gives

    <W_mu(t)> <=C_mu exp(v_mu|t|), v_mu=12(exp(mu)-1).

The complete tail changes still have norm <=g_max=1.1. Assisted prefix pop,
bare matter-driven end motion and cut-crossing interactions are the only
uncanceled drives. Conversion sources lie inside b; no source is hidden
in the empty tail. With L=ell-2-b=3N/16-1, mu=1,

    ||(U_full(t)-U_prefix(t))psi||<=D_A(t),
    D_A(t)=min{2,tau_(a+1)(1.1|t|)
           +(14 sqrt(C_1)/v_1)exp(-L/2)[exp(v_1|t|/2)-1]}.    (7)

This is rederived using resource q_x, not inherited by renaming P1's
number-conserving operating theorem. All full-Hamiltonian excursions and
returns occur on the physical side of exact Duhamel. At N>=2^18 and
0<=t<=pi/2+N/65536, elementary e bounds give

    D_A(t)<=exp(-N/5)+exp(-N/12).                            (8)

For example v_1<24, C_1<8, prefactor<7/2 and t<2+N/65536. The second
exponent is at most 24.5-(3/32-1/8192)N, which absorbs its prefactor into
exp(-N/12) already for N>=4096. The head tail is smaller than exp(-N/5).

Any event probability differs by at most 2D_A. Hence (6) becomes

    Pr_full(P(pi/2)>=K)>0.946-2D_A(pi/2),
    Pr_full(P(pi)>=K)<0.215+2D_A(pi)<1/4.                    (9)

The frozen use window includes pi. A single time violating its universal
requirement refutes that gate for the entire class containing this point.
It does not prove permanent absence of matter, eventual thermal behavior,
or failure of every future fuel/bath/clock architecture. In particular do
not reinterpret a time-averaged or later-window statement as a pass of this
every-time, already frozen protocol.

## 4. Locality and energy survive the failed retention gate

The same J_0=6 connected interaction chains give a prefix information bound

    ||[S,U_prefix^*OU_prefix]||
       <=2||S||||O||min(1,s exp(-mu d)[exp(24e^mu|t|)-1]),

for disjoint local physical supports, with s the smaller support size.
On-site mass/fuel energies rotate operators locally; they are not globally
ignored. Full-state comparison adds 4||S||||O||D_A to source response.
Prefix-loss and storage-error probabilities are <=D_A^2 when initially
consistent. Thus the counterexample is not an information shortcut or an
energy blowup: an operationally local, extensive-energy medium can still
return matter to its prepared fuel. Initial O(N) resources and local motion
alone are insufficient to force a sustained matter-bearing regime.

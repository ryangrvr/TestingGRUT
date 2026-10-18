# Exact operational reduction of the successful instrument

**Tier:** author-derived algebraic theorem; independent review pending.
It does not assert a small moving-geometry error or operating pass.

## 1. The experiment is unchanged

Use the B1 moving law and the B2-RB conditional preparation: a given
recorded word w of length L in 24N, L>=2016, q=a_w, actual conditional
64N^2-mode receiver sigma_w with L-1 bosons, vacuum charges/fields and
excited fuels at the recorded root and right label b=w_(L-1). All inactive
fuel slots remain in their exactly invariant vacuum. The two inputs are
arbitrary qubits, possibly entangled with reference Z, independent of the
spatial preparation. Fuel conversion/recombination and geometry stay on.
Retain the original accepted-record restrictions on N,L and the obtained
sigma_w; no additional record selection is introduced. These algebraic
identities do not themselves invoke or validate the formation theorem.

R is the SAME final success projection as B2-RB: the original recorded
prefix/bonds/storage intact, both fuels empty, co-pairs at one recorded
site in X_L=[L/8,3L/8-1] and one in Y_L=[5L/8,7L/8-1]. Decode plus
colors in recorded left/right order; trace spatial, tail and receiver
registers. Every other outcome is erasure. Grown words may satisfy R.
No geometry history is measured or selected.

T is the SAME retention event: original prefix/storage survival and k=2,
including separated constituents. Output uses I_L=[L/(8J),L/(6J)];
retention uses [2/J,L/(6J)], J=1/3840. Thresholds remain 1/10 output,
1/4 retention and 1/100 conditional trace distance. There is no new
probability claim about selecting a suitable record or preparing fuel.

## 2. Complete Gram description

Encode the two color carriers in current word order, with root fuel before
its plus register and right fuel after its plus register. The exact
color-blind Hamiltonian has H=A_c tensor I+B_c tensor SWAP. These are
symmetry blocks of ONE law, not alternative physical Hamiltonians. Put
H_+=A_c+B_c, H_-=A_c-B_c, U_s=exp(-itH_s), and purify the ACTUAL sigma_w
in the initial spatial/environment vector xi_w. Let

    chi_s=R U_s xi_w,
    u=||chi_+||^2, v=||chi_-||^2, c=<chi_-,chi_+>.

Then 0<=u,v<=1 and |c|^2<=uv. For Pi_+,Pi_- the color triplet/singlet
projections, the full UNNORMALIZED successful operation is

    E_R(rho)=u Pi_+rho Pi_+ +v Pi_-rho Pi_-
                   +c Pi_+rho Pi_- +conj(c) Pi_-rho Pi_+.        (1)

All formulas below also act on arbitrary external references. No successful
probability has been divided out. In the fixed-path experiment the same
operation is p_f(t) rho, where p_f is the ACTUAL fixed-path probability,
not an adjustable value or merely its lower bound.

For completeness, after erasing the unread output qubits the FAILURE
reference state is

    (1-u)rho_Z,+ +(1-v)rho_Z,-,
    rho_Z,s=Tr_AB[(Pi_s tensor I_Z)rho(Pi_s tensor I_Z)].

It is not generally (1-p_R(rho))rho_Z: input-dependent success can steer
the reference on either flag. The fixed-path u=v=p_f case does factor.

For a pure input including a reference, set x=||Pi_+psi||^2. Its two
orthogonal components reduce (1) to a 2-by-2 Hermitian matrix on their
span. Every x in [0,1] is attained by a triplet/singlet superposition
without a reference. Purification and partial trace extend the resulting
bounds to all mixed/reference inputs. The full input quantifier is thus
solved, without assuming a classical mixture of parity sectors.

## 3. Exact unnormalized operation distance

Define eta_R=sup_(rho_ABZ) ||[(E_R-p_f id) tensor id_Z](rho_ABZ)||_1,
with NO factor of one half. Let a=u-p_f, d=v-p_f, z=c-p_f. On a pure
input the difference has matrix

    [ a x,                  z sqrt(x(1-x)) ]
    [ conj(z)sqrt(x(1-x)),   d(1-x)        ].

Its trace norm is the maximum of absolute trace and eigenvalue gap.
Since the trace maximum is already attained at an endpoint, the EXACT
reference-uniform distance is

    eta_R^2 = max_(0<=x<=1)
                  [(a x-d(1-x))^2+4|z|^2 x(1-x)].              (2)

This is a single quadratic. Its endpoints give d^2,a^2. If
K=(a+d)^2-4|z|^2<0, also test the vertex

    x_*=[d(a+d)-2|z|^2]/K

when it lies in (0,1); otherwise endpoints suffice. References cannot
increase (2), and a pure color input attains its maximum. The theorem
requires no unknown dimension-dependent constant or general-state norm.

If the entire flagged success/erasure instrument is compared, its exact
reference-uniform state trace norm is

    max_x [N(x)+|a|x+|d|(1-x)],

where N(x) is the 2-by-2 matrix trace norm above. Orthogonal reference
flags on the two parity components attain the failure contribution. This
is distinct from eta_R, which deliberately measures the successful map.

## 4. Exact output and conditional information criteria

The success probability on the same input is q(x)=u x+v(1-x). Hence

    p_R,min=min(u,v).                                          (3)

For u,v>0, the EXACT worst conditional trace distance to the original
input, with trace distance defined as half the trace norm, is

    D_R^2=max_(0<=x<=1)
       x(1-x)[(1-x)|c-v|^2+x|c-u|^2]/[u x+v(1-x)]^2.           (4)

For mixed inputs, purify first, condition on the same success probability
and trace the extra reference; thus the same exact supremum applies.
No probability or channel quality is inferred from the other.

Equation (4) turns the 1% input-uniform condition into positivity of a
specific cubic on [0,1]. Set e0=|c-v|^2, e1=|c-u|^2, d0=u-v and
epsilon=1/100. The condition is NECESSARY AND SUFFICIENT:

    F(x)=epsilon^2[v+d0 x]^2
              -x(1-x)[e0+(e1-e0)x] >=0 for every x in [0,1].    (5)

Its coefficients are epsilon^2 v^2, 2epsilon^2 v d0-e0,
epsilon^2 d0^2-e1+2e0, e1-e0. Given controlled Gram entries, only endpoints
and roots of a quadratic derivative need be checked. This finite input
test does NOT solve the remaining time/size/receiver quantifiers.

For intuition, if u=v=p>0, D_R=|p-c|/(2p). Perfectly equal success
probabilities can coexist with conditional error 1/2 (c=0) or 1 (c=-p).
These are channel examples, not counterexamples to the actual B1 law.

An additional useful identity is, with e^2=u+v-2Re(c) and D=uv-|c|^2,

    D_R(x)^2=x(1-x)[e^2 q(x)-D]/q(x)^2
                <= e^2/[4 min(u,v)].                          (6)

Thus a bound on outgoing-vector difference can control information
without preserving common geometry. INFLUENCE_REDUCTION.md identifies
the actual restricted coupling that generates this difference.

## 5. Retention is separate and simpler

Let u_T=||T U_+xi_w||^2 and v_T=||T U_-xi_w||^2. The retention probability
is u_T x+v_T(1-x); its worst value is min(u_T,v_T). Cross-parity coherence
does not enter this scalar probability. Compared with the actual fixed
retention p_f,T, its reference-uniform probability error is exactly

    eta_T=max(|u_T-p_f,T|,|v_T-p_f,T|).                         (7)

R is a subevent of T, so u<=u_T and v<=v_T, but the different 10%/25%
requirements prevent replacing retention by output success.

## 6. Two honest ways to close the physical gate

Comparison route: B2-RA gives p_f>p0=514089/4000000 and
p_f,T>q0=434281/1000000 over their respective windows. If uniform
eta_R,eta_T are proved, then

    p_R,min>p0-eta_R, p_T,min>q0-eta_T,
    D_R<=eta_R/(p0-eta_R) when eta_R<p0.                       (8)

The last bound follows by normalizing only AFTER bounding the successful
map: trace error <=eta_R and normalization costs at most another eta_R.
The factor one half in trace distance leaves the ratio in (8).
The sufficient success-map budget for 1% is eta_R<=p0/101. At
eta_R=eta_T=1/1000, this would give output >0.12752225, retention
>0.433281 and D_R<=4000/510089<1/100. These are UNACHIEVED certificate
values. They are operation-norm errors, not B2-RB's state-vector errors.

Direct route: no closeness to p_f is needed if the ACTUAL moving entries
obey min(u,v)>=1/10, min(u_T,v_T)>=1/4 and (5) throughout the declared
windows. These are necessary and sufficient operational conditions for
the unchanged task. A larger common success probability need not resemble
the fixed probability to be useful. A rare perfect channel still fails
the probability gate.

Neither route is completed here: the moving entries and their continuous
window bounds remain unknown. The result is an exact operational reduction,
not another scalar source surrogate or an operating pass.

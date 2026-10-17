# Which receiver assumptions do the formation work?

This follows the S9 rederivation; it is not independent validation of it.
The benchmark law, L=64N^2 and its failed numerical branch are unchanged.
The lemma below is a same-author generalization. Its example families are
separately declared receiver variants, not observations of unchanged WM2.

## A sufficient correlation criterion

Keep the complete word algebra, lambda=g/N, initial receiver vacuum and
g in a fixed positive compact interval. Fix the original time protocol
4N^2<=t<=8N^2. For a finite self-adjoint receiver matrix h_N, suppose
its contact correlation C_N admits a reference correlation C_* with
constants independent of N, disorder and any receiver-family parameter:

1. |C_*(u)|<=A(1+u)^(-p), p>1, for u>=0.
2. Gamma(delta)=2 Re integral_0^infty C_*(u)exp(-i delta u)du
   is bounded below by gamma>0 over the whole detuning box.
3. sup_(0<=u<=8N^2)|C_N(u)-C_*(u)|<=epsilon_N and
   N^2 epsilon_N tends to zero.

Then the current proof remains history-retaining. Its contact bound is

    b_N=lambda sqrt(M)(integral |C_*|+8N^2 epsilon_N)=O(N^(-1/2)),

and the two-commutator bound B_N is O(1). Write

    R_N=integral_0^(8N^2) |C_*(u)| min(B_N u,2sqrt(M)) du,
    L_*(t)=integral_t^infty |C_*(u)|du.

The exact-current remainder is bounded by

    2lambda^2 sqrt(M) R_N
      +2lambda^2 M L_*(t)+2lambda^2 M (8N^2)epsilon_N.             (R1)

This expression follows by the same exact adjacent-charge identity and
retardation derivative, with no assumption about commuting histories or
independent emissions. Splitting R_N at sqrt(N) gives

| Kernel exponent | R_N | Contribution to f/N after damping |
|---|---|---|
| 1<p<2 | O(N^((2-p)/2)) | O(N^((1-p)/2)) |
| p=2 | O(log N) | O(N^(-1/2) log N) |
| p>2 | O(1) | O(N^(-1/2)) |

The normalized missing-tail contribution after integration is
O(N^(2-2p)) for 1<p<2, O(N^(-2) log N) for p=2 and O(N^(-2)) for
p>2. The normalized finite-correlation contribution is O(N^2 epsilon_N).
All vanish. Gronwall and the same free-count event inequality therefore
give, uniformly over the window,

    1-p_N(t)<=2 exp(-4 gamma g_min^2)+o(1).                       (R2)

In particular if 2exp(-4 gamma g_min^2)<1/2, macroscopic formation is
>1/2 for all sufficiently large N. If gamma is smaller, this criterion
alone gives a weaker or vacuous result for the frozen time protocol;
one may not retime it afterward and call the original test passed.

R1-R2 broaden sufficient receiver properties, not the meaning of a phase.
The explicit finite correlation condition and its uniform constants are
obligations, not consequences of saying a receiver has many modes.

## A nonshrinking structured family passes this sufficient criterion

Define the distinct receiver family R by

    h_N(a,omega)=omega I+a h_path,L,
    L=64N^2, a in [9/10,11/10], |omega|<=1/4.

All other coefficients, preparation and the observation interval stay
fixed. The positive physical frequencies 4+sigma(h_N) remain positive.
The limiting kernel is exp(-i omega u) C_infty(a u), and

    Gamma(delta)=sqrt(4a^2-(delta+omega)^2)/a^2.

For |delta|<=3/4 this is >=sqrt(4(9/10)^2-1)/(11/10)^2>6/5.
Its decay exponent is p=3/2, uniformly in a,omega. A uniform finite
closed-walk bound is

    |C_N-C_*|<=(160/69)(33/80)^(128N^2)<2^(-120N^2).

The scalar omega phase cancels in the finite/infinite difference. The
ratio of successive Taylor terms is <=11/80 and e<3 supplies 33/80.
Using gamma=6/5, g_min=9/10 in R2 gives a limiting failure upper bound
2exp(-486/125)<1/2. Thus this family also forms macroscopically at
every time in the same window for sufficiently large N, uniformly in
its receiver parameters and the storage disorder. This does not assert
S9's same finite threshold 2^28 for the family; no new threshold was
optimized or certified here. It does not reduce mode-count scaling.

## An equal-mode-count receiver can fail at every time

Define a separate flat receiver control h_flat=2I on the same L modes.
Keep the same word, storage, coupling and initial state. It has positive
physical frequencies 6, and no spatial adjacency. In the reachable
vacuum-seed sector all emitted photons remain at contact and n_R=K.
Its diagonal satisfies

    D_flat=sum_v (2+delta_v) P_(v attached)>=(5/4)K.

The assisted Hermitian interaction still has norm <=g. The seed's mean
energy is zero and conserved, so for every time

    (5/4)<K><=<D_flat>=-<V><=g,
    Pr[m>=ceil(N/2)]<=4g/[5(ceil(N/2)-1)], N>=4.                (R3)

The bound tends to zero. The large receiver count, vacuum, reversibility,
positive frequencies and label symmetry alone do not ensure formation.
This is a resource-dependence counterexample in a changed receiver law,
not a failure inside S9's frozen path receiver.

## Preparation, normalization and resource conclusions

The declared bondless consistent seed has q=root only. The exact charge
n_R=|q|-1 therefore forces receiver vacuum: it is not an independently
adjustable state within that seed. Changing it requires a new complete
preparation and energy ledger. With a nonvacuum preparation the homogeneous
field term in F1 generally survives; this proof supplies no robustness
theorem for it.

The balance lambda=g/N, t_form=O(N^2) gives nonvanishing gamma lambda^2 t.
Changing normalization requires rechecking both that time scale and the
commutator remainder. A longer run is not the unchanged acceptance protocol.
For the path receiver, L proportional to the O(N^2) horizon gives a large
no-return buffer in the current walk estimate. It is sufficient; no lower
bound proving L=Omega(N^2) is physically necessary has been established.

The construction tolerates a finite family of bandwidth/center changes,
but does not explain its receiver, bosonic statistics, charge capacity or
vacuum-seed convention. Lowering a sufficient N threshold and reducing
these supplied physical resources remain different questions.

# Charge-balanced formation-and-preparation estimate

Author-derived; independent review pending. CHARTER.md fixes the protocol.
S9 is a dependency for the original-trial corollary only. The conditional
preparation proof does not require S9's formation estimate.

## 1. Exact local conversion, including reverse motion

On an intact selected bond, with no operating hopping during the pulses,
the relevant local states are |0> (empty fuel/matter/field),
|F> (one fuel, no matter, E=0), and |M> (fuel empty, + at parent, - at child,
E=+1). Their diagonal energies are 0, Omega=3+Delta, and 3 respectively.
The last energy is two mass units plus one electric unit.

The loading rotating frame gives Delta |F><F|+eta(|F><0|+adjoint).
The fixed load time pi/2 produces

    p_L(Delta,eta)=eta^2/(eta^2+Delta^2/4)
                    sin^2[(pi/2)sqrt(eta^2+Delta^2/4)].       (6)

Conversion is the exact Hermitian two-level matrix on {|F>,|M>},

    [[3+Delta,lambda],[lambda,3]],

with |0> dark, giving the same expression p_C(Delta,lambda).
Final matter-pair probability is p_e=p_L p_C. Relative phases and fuel
coherences are retained, although number-event probabilities do not depend
on those phases. Recombination is part of this solution, not a discarded
correction.

For |Delta|<=1/2 and a in [0.9,1.1],

    a^2/(a^2+Delta^2/4)>=324/349,
    0.9<=sqrt(a^2+Delta^2/4)<8/7.

The sine's angle lies within pi/14 of pi/2. Using sin^2(pi/2+x)>=1-x^2
and pi<22/7 gives

    p_L,p_C >=(324/349)(2280/2401)>7/8,
    p_e >=49/64.                                             (7)

This works for independently different parameters on every selected bond;
no probability distribution or parameter-dependent pulse retiming is needed.

## 2. Common-space comparison with moving geometry

Choose ell=N/4 and retain its recorded prefix. Define a self-adjoint reference
on the full physical Hilbert space: on the subspace retaining that prefix,
reflect assisted geometry at length ell, omit bare geometry there, and retain
the entire receiver, storage diagonal, drive and conversion. On its orthogonal
complement use the compression of the actual stage Hamiltonian, with reference
cross blocks zero. The physical Hamiltonian itself is unchanged.

In this reference the matching lies wholly inside the retained prefix and
all tail matter/fields/fuel are vacuum. Prefix storage is a spectator. Thus
the K independent three-state pulse systems factor from tail/receiver evolution,
even for the actual conditional receiver state. Reflection is a mathematical
comparison, not a physical clamp. The initial word has m>=N/2, and
a=m-ell>=N/4 assisted moves are needed to reach the reflecting end.

The whole changing tail operator has norm <=g_max=11/10. All extensive
diagonals, prefix pulses and receiver-only motion preserve tail length.
The complete Dyson series therefore bounds its end amplitude by
tau_a(g_max t), tau_r(x)=sum_(n>=r)x^n/n!. The only noncanceling drive
between actual and reference during these pulses is the assisted prefix-
breaking pop. Bare terms vanish on the reference: its head/parent lies
in the empty region, since 2K=N/16<ell-1. Exact Duhamel over total load
plus conversion time pi gives the **full-state** bound

    ||psi_prepared-psi_reference||<=e_p,
    e_p=tau_(a+1)(11pi/10) <=exp(-N/4), N>=256.               (8)

For the coarse inequality use 11pi/10<4, k=a+1>N/4, and
tau_k(4)<=e^4(4e/k)^k. At k>=65 this is <=exp(-N/4); the ratio improves
with k. This charges all physical geometry/storage/receiver excursions
and their returns, not just the final omitted amplitude. No electric or
matter configuration is removed from the physical Hamiltonian.

Loading itself is exactly independent of geometry even if the prefix moves:
its drive acts on the chosen potential-pair fuel registers, and no matter
exists yet. Equation (8) is a conservative unified bound for both intervals.

## 3. Extensive matter probability with the buffer paid

The reference number event is a sum X of K independent Bernoulli variables
with success probabilities p_e>=49/64, by the exact product pulse state.
Independence here is an exact tensor-factor fact for the comparison during
preparation, not an approximation to interacting operation or receiver emission.
Hoeffding's elementary bounded-variable estimate gives

    Pr_ref(X<K/2)<=exp[-2K(49/64-1/2)^2]
                   =exp(-289K/2048).                        (9)

The same conclusion follows by the exponential-moment product bound.
The non-strict endpoint K/2 only makes this bound more conservative.
There are P=2X particles; K is even. A projector probability changes by
at most 2e_p for normalized pure states separated by e_p (purify mixtures).
Consequently the **actual** preparation obeys, uniformly over the frozen
parameter class and conditional receiver state,

    Pr(P>=N/32 | accepted record)
      >=1-exp(-289N/65536)-2exp(-N/4).                        (10)

There is no particle-number postselection in the protocol. Unconverted fuel
and failed conversion outcomes remain in the state. All reference matter is
supported at b=2K-1=N/16-1 by construction. The empty region through ell=N/4
is therefore produced by the same paid controls, not obtained by a separate
favorable buffer measurement. Physical deviations from that support are
charged in e_p. Once operating begins, N_+,N_- are conserved, so (10) persists.

At N=256 this elementary conditional lower bound is already greater than
0.676. This is a formula evaluation, not an N=256 simulation or formation
guarantee. No useful S9 formation statement at that N is implied.

## 4. Original trials and exactly what is inherited

Before recording, all fuel and matter are vacuum and every preparation
control is off. The new positive field energy and empty-head gates have zero
effect there. The original WM2 formation Hamiltonian is exactly recovered.
If the archived S9 theorem is valid, its uniform statement at N>=2^28 and
T_f=4N^2 gives p_form>1/2, with original root weights, L_N=64N^2 and the
declared independent storage disorder. Multiplying (10),

    Pr(accepted AND P>=N/32)
      >(1/2)[1-exp(-289N/65536)-2exp(-N/4)]>1/4.             (11)

N must also be a multiple of 64. Rejected formation trials count zero. The
recording intervention and work source are supplied; controller reliability
and laboratory implementation are not inferred from this theorem. No receiver
reset is used. At the inherited threshold the receiver still has 2^62 modes
and formation time is 2^58. P1 does not reduce those resources.

Equations (6)-(10) and the conditional use theorem are new author derivations.
Equation (11) additionally depends on S9, itself author-derived and pending
independent proof review. It must not be described as independently validated.

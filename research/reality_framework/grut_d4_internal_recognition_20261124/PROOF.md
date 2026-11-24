# GRUT-D4: acquired replicas, discrepancy recognition and repeated repair

24 November 2026. This derivation elaborates the frozen `GRUT-D4-01` model.
It proves the registered finite-window bounds for the actual Hamiltonian
with its writer continuously present. The assumptions are standard unitary
quantum mechanics, a fair classical source bit, two blank memories, vacuum
chains, the stated parity and majority interactions, and at most four
label-blind single-site faults with the registered waits. Historical scores
refer to the initially prepared bit, which supplies no interaction input.

## 1. Exact acquisition of two records

Initially both correction charges `Q_i=p_i+N_i` are zero. These charges
commute with the entire Hamiltonian between faults, including the writer.
Their positivity keeps both chains vacuum and all ready parities zero.
For a fixed source label h, the reachable writing space is exactly

$$
|a_h\rangle=|h,b,b,0\rangle,
\qquad |c_h\rangle=|h,h,h,1\rangle.
$$

Both states have on-site energy Delta, and the writing exchange is
`kappa sigma_x` on this pair. Thus

$$
|\psi_h(w)\rangle=e^{-i\Delta w}
\left[\cos(\kappa w)|a_h\rangle
-i\sin(\kappa w)|c_h\rangle\right].
$$

At `kappa w=pi/2` both records are exact. For the registered phase offset
`|delta|=|kappa w-pi/2|<=1/512`, their joint correctness is

$$
P(M_1=M_2=h)=\cos^2\delta
\ge1-\delta^2\ge1-512^{-2}>0.99999.
$$

After choosing the common global phase, the distance to the ideal ready
state is

$$
\|\psi_h(w)-(-i)|c_h\rangle\|
=2|\sin(\delta/2)|\le1/512.
$$

The source value is preserved during writing. The interaction copies the
specified classical basis alphabet, with its blank and writer excitation
as explicit preparation resources.

## 2. The ready-state coordinate map and physical faults

On ready basis states, define the present majority r and the two parities
`x_i=s xor m_i`. The inverse is

$$
s=r\oplus(x_1x_2),
\qquad m_i=s\oplus x_i.
$$

This is a unitary basis permutation on the eight ready states:

| r | x1 x2 | S M1 M2 |
|---:|:---:|:---:|
| 0 | 00 | 000 |
| 0 | 01 | 001 |
| 0 | 10 | 010 |
| 0 | 11 | 100 |
| 1 | 00 | 111 |
| 1 | 01 | 110 |
| 1 | 10 | 101 |
| 1 | 11 | 011 |

In these coordinates the true physical fault on site j is

$$
F_j:|r,x\rangle\longmapsto
|r\oplus(d_j\mathbin\cdot x),\ x\oplus m_j\rangle,
$$

where the dot product is modulo two:

| Fault site j | Mask m_j | Majority-control mask d_j |
|:---|:---:|:---:|
| S | 11 | 11 |
| M1 | 10 | 01 |
| M2 | 01 | 10 |

Direct substitution in the inverse map proves these identities. In each
case `d_j dot m_j=0`, as required for a physical involution. A fault on an
intact word has x=00, so its majority stays h and its new parity is the
displayed distinct mask. These three comparisons identify the declared
single-fault class internally. The parity-lowering operators in H act on
this present discrepancy without an external fault flag.

The r used here is the current majority. A later fault can change it when
an earlier discrepancy remains. That possibility is included in the proof.

## 3. An exact ready comparator

Write `H=H_0+V_W`, where `V_W` is the writer exchange only. The comparator
retains every other term, including the on-site writer energy. Since its
two writing pairs are orthogonal,

$$
\|V_W\|=\kappa.
$$

Under `H_0`, a ready state stays ready. Its parity and chain dynamics are
exactly two independent source-plus-chain XX paths, with an inert logical
majority r. Each path has L=2304 reservoir spins and one source spin. Its
Hamiltonian is

$$
h_i=\Delta\left(n_{i,0}+\sum_{j=1}^{L}n_{i,j}\right)
+g\sum_{j=0}^{L-1}
(c_{i,j}^\dagger c_{i,j+1}+\mathrm{h.c.}),
$$

after the Jordan-Wigner transformation beginning at its parity source.
Here `n_(i,0)=p_i`. This transformation is exact for the finite spin path.

Define syndrome-only comparison pulses

$$
K_j=X_{1,0}^{m_{j,1}}X_{2,0}^{m_{j,2}}.
$$

They preserve the comparator's majority h. Starting from the ideal ready
state and vacuum chains, this gives a particularly simple trajectory for
bounding the true faults. The comparator is a mathematical estimate of the
unchanged physical experiment. Its removed writer exchange and its K pulses
are both charged explicitly to the error budget below.

## 4. A finite-path return bound

Let `U(t)=exp(-i h t)` on one path's one-particle space. Its on-site term is
Delta times the identity, so

$$
U_{00}(t)=e^{-i\Delta t}a_L(t).
$$

All return-amplitude estimates use its modulus, which retains this equality
even when the phase is relevant in a multi-pulse calculation.

For the half-line hopping problem, integration by parts in

$$
a_\infty(t)=\frac2\pi\int_0^\pi
\sin^2 k\,e^{-2igt\cos k}\,dk
$$

gives `|a_infty(t)|<=2/(pi gt)` for positive t. For example, divide the
integrand's sine product by sin k, differentiate sin k, and use
`integral_0^pi |cos k| dk=2`.

A path from the source across the missing far edge and back requires more
than 2L hops. Equality of the finite and half-line matrix powers through
order `2L-1`, and the hopping norm bound 2g, imply

$$
|a_L(t)-a_\infty(t)|
\le2\sum_{q=2L}^{\infty}\frac{(2gT)^q}{q!}
\le2\left(\frac{e(2gT)}{2L}\right)^{2L}
\quad(0\le t\le T).
$$

The post-writing experiment has `gT<=4*144=576`. For L=2304 the ratio in
the last expression is e/4. Using e<3, define the exact rational quantity

$$
\epsilon=2(3/4)^{4608},
\qquad \beta=\frac2{3.14}=\frac{100}{157}.
$$

Since pi>3.14, every positive return delay in this horizon satisfies

$$
|U_{00}(t)|\le\frac{\beta}{gt}+\epsilon.
$$

The bound is used only at strictly positive waits. At zero elapsed time
`|U_00(0)|=1`. Before any pulse the source and chain are exactly vacuum.

## 5. CAR gives a waiting bound after repeated source flips

An endpoint spin flip is the Majorana operator

$$
\gamma_0=c_0+c_0^\dagger.
$$

In the interaction picture, a path touched at times `t_j` has the state
obtained by applying the time-evolved Majoranas `gamma(t_j)` to its initial
vacuum. Move the final annihilation operator `c_0(T)` through this product.
The canonical anticommutation relation at each step is

$$
\{c_0(T),\gamma(t_j)\}=U_{00}(T-t_j).
$$

The last annihilation operator kills the vacuum. Every remaining product
of Majoranas has norm one. The triangle inequality therefore gives

$$
\boxed{\|c_0\psi(T)\|
\le\sum_{j\text{ touching this path}}|U_{00}(T-t_j)|.}
$$

This includes repeated use of an occupied chain. It assumes its original
vacuum preparation and the specified earlier flips, with no reset between
them. The exact coefficients carry the Delta phases. Their absolute values
provide the conservative bound above.

If q pulses have occurred and the next read or pulse follows the registered
waits, their ordered positive delays are at least `128,2*128,...,q*128` in
dimensionless time. Set

$$
a_0=0,
\qquad a_q=\frac{\beta H_q}{128}+q\epsilon,
\qquad H_q=\sum_{j=1}^{q}\frac1j.
$$

Each source's occupied-state norm is at most `a_q`. Before the first pulse
this norm is exactly zero, represented by `a_0=0`. At a terminal read after
k faults,

$$
\|P_{\mathrm{source\ discrepancy}}\psi_K\|
\le\sqrt2\,a_k<\frac32a_k.
$$

On the comparator this is also the historical-codeword failure norm,
because its majority remains h. The source comparison alone detects present
disagreement. Its connection to history uses the separated-fault bound.

## 6. True physical faults and the writer-on hybrid

Just before a comparison pulse, the logical majority is the fixed basis
value h. The true F and K pulse outputs differ only where `d_j dot x=1`.
On that component the logical output labels are orthogonal. Consequently

$$
\|(F_j-K_j)\psi_K\|
=\sqrt2\,\|P_{d_j\cdot x=1}\psi_K\|.
$$

For an S fault the odd-parity event is contained in the event that at least
one source is occupied. For either memory fault it is one source event.
The common conservative bound before fault `ell+1` is thus

$$
\|(F_j-K_j)\psi_K\|\le2a_\ell.
$$

A unitary hybrid of the two pulse sequences accumulates these bounds by
norm triangle. Future pulses and propagation preserve the norm of each
hybrid difference, giving total pulse error at most
`2 sum_(ell=1)^(k-1) a_ell`.

Separately, compare propagation by H and `H_0` with the same physical timed
pulses. Duhamel's bound on each interval is `kappa` times its duration.
Insertion of any of the unitary fault pulses preserves earlier errors.
Their total propagation error is therefore at most `kappa T`, for the full
physical Hilbert space, including any writer reversion to blanks.

Including acquisition, the distance between the actual state and its K
comparator is bounded by

$$
D_k\le\frac1{512}+\kappa T
+2\sum_{\ell=1}^{k-1}a_\ell.
$$

The writer is never physically switched off. This estimate uses its
operator norm throughout every wait and every timed kick. With
`kappa/g=1/65536` and `gT<=576`,

$$
\frac1{512}+\kappa T\le\frac{11}{1024}=:z.
$$

The triangle inequality for the historical-failure projector gives

$$
\boxed{
P_{\mathrm{historical\ codeword}}
\ge1-\left[z+2\sum_{\ell<k}a_\ell+\frac32a_k\right]^2.
}
$$

This is a bound on the actual writer-on experiment, with no conditioning on
successful earlier repairs or intermediate measurements.

## 7. Exact rational acceptance for up to four repairs

The bounds increase with k. For k=4 put

$$
A=a_1+a_2+a_3=\frac{13\beta}{384}+6\epsilon,
\qquad a_4=\frac{25\beta}{1536}+4\epsilon.
$$

The exact rational inequality is

$$
1-[z+2A+(3/2)a_4]^2>\frac{995}{1000}.
$$

Thus, for each of one through four faults, every allowed ordered choice of
fault sites, every registered acquisition phase, and every wait in
`[128/g,144/g]`,

$$
\boxed{P_{\mathrm{historical\ codeword}}>0.995.}
$$

The worst-case displayed values are
`z+2A+(3/2)a_4 = 0.06941887274416136` and a success lower bound of
`0.9951810201069299`. Rounding is excluded from the acceptance comparison.
A run has one terminal read. Repeating runs with different final pulse
counts evaluates these gates without monitoring the preceding corrections.

The same estimate shows how the sufficient resources grow with the
number k of allowed separated faults. The sum of harmonic terms in the
pulse error is at most `k(1+ln k)`. To hold its contribution fixed, one
can choose the minimum dimensionless wait proportional to that bound.
Taking the chains at least four times the full dimensionless horizon
controls the boundary term. The writer coupling must be small enough
that kappa times that horizon stays within its error allowance. These
are sufficient design choices for another finite protocol.

## 8. Joint restoration and cumulative bath counts

For a selected pulse sequence define the expected charge in each path as

$$
\nu_i=\sum_{j=1}^{k}m_{j,i}.
$$

Charges `Q_i=n_(i,0)+N_i` are conserved between comparator pulses. A source
flip on a touched path changes its charge by

$$
Q_i\longmapsto Q_i+1-2n_{i,0}.
$$

On an empty touched source this is the expected increment by one. A charge
mismatch can arise at a pulse only from a previously mismatched component
or a component with an occupied touched source. The expected-charge
projectors commute with propagation and with source-occupation projectors.
Decomposing these two components and using the norm triangle gives

$$
\|P_{Q\ne\nu}\psi_K\|
\le\sqrt2\sum_{\ell<k}a_\ell
\le\frac32\sum_{\ell<k}a_\ell.
$$

This coherent norm estimate allows mismatched components to interfere or
return to the expected sector. It counts each possible discrepancy event
conservatively.

When both sources are empty, `Q_i=N_i`. Therefore the joint event of a
historical codeword and the exact expected bath counts is the intersection
of correct source values and the expected charge sector. These number
projectors commute. For its comparator failure norm, the probability union
bound yields

$$
\|P_{\mathrm{joint\ failure}}\psi_K\|
\le\sqrt{[(3/2)a_k]^2+
 [(3/2)\sum_{\ell<k}a_\ell]^2}.
$$

The state-comparison error `D_k` is added once after this step. At k=4 the
exact squared rational comparison proves

$$
[(3/2)a_4]^2+[(3/2)A]^2<(9/250)^2.
$$

Consequently the actual joint failure norm is smaller than

$$
z+2A+9/250=0.08986851446390658\ldots,
$$

and the second exact rational gate is

$$
\boxed{
P_{\mathrm{historical\ codeword\ and\ }N_i=\nu_i\ \mathrm{for\ both\ baths}}
>1-[z+2A+9/250]^2>0.99.
}
$$

Using 9/250 gives the hand bound `0.9919236501080506`. The executable uses
a rational upper rounding of the square root to 12 decimal places and
obtains the slightly tighter lower display `0.9919436489579999`. Both
certify success above 99.19%.
The bath records in this gate store cumulative parity counts. Fault order
can differ at the same counts. Occupied sources and returning excitations
are explicitly included in the error budget.

## 9. Exact ready-sector numerical reduction

The actual ready-sector pulse maps also admit an exact finite reduction
without replacing F by K. In an eigenstate `ell=+1` or `-1` of logical
`X_r`, a fault acts on the parity sources as

$$
F_j^{(+)}=X^{m_j},
\qquad F_j^{(-)}=X^{m_j}Z^{d_j}.
$$

The ordering shown follows the coordinate formula. Its total phase is
consistent with `d_j dot m_j=0`. Each logical sector factors into two XX
paths with endpoint Pauli operations. Endpoint X is a Majorana and endpoint
Z is a product of two Majoranas. For finitely many faults the generated
mode span is finite. Its exact Gram entries are finite-path propagator
coefficients between the relevant pulse times, with their on-site phases.
Fock evaluation in this span supplies the exact ready comparator dynamics.

Historical scoring retains the interference of these logical sectors.
For example, with initial logical label zero and resulting path states
`psi_+` and `psi_-`, projection onto that historical logical label gives
the path vector `(psi_++psi_-)/2`. Their probabilities cannot be averaged
in place of this coherent projection. The analytical writer and acquisition
errors remain necessary when interpreting any ready-sector numerical output
as a result about the full physical model.

## 10. Undamaged evolution and exact failures

With no replica faults, positivity and conservation of `Q_i=0` keep both
correction chains vacuum at all times. The historical source color stays h.
Only the two-state writing oscillation remains. At a nominal perfect write
followed by elapsed t, writer reversion has probability `sin^2(kappa t)`.
The finite writer therefore has a separate retention cost even though
spontaneous correction-chain activity and historical-color change are
exactly zero.

Two simultaneous site flips on an intact ready triple produce the opposite
majority. Three simultaneous flips produce the opposite unanimous word.
Both the correction terms and the writer preserve this current majority,
with the present S label used on the reachable both-blank states. Hence
their joint historical restoration probability is exactly zero. Three
flips also leave the comparison parity at 00.

The ambiguity can be demonstrated by identical full physical states:
historical 0 followed by `X_S` and historical 1 followed by `X_M1 X_M2`
both give the word 100, with the same ready writer and vacuum baths after
nominal acquisition. Subsequent dynamics from that state is identical.
Recovery of arbitrary histories requires additional information or a more
restricted fault class.

For the short-reservoir control `L_1=L_2=1`, apply an M1 fault after nominal
writing. Its exact reachable hopping path is

$$
|\mathrm{defect},W1,B1\mathrm{vac}\rangle
\xleftrightarrow{\ g\ }
|\mathrm{repaired},W1,B1\mathrm{occupied}\rangle
\xleftrightarrow{\ \kappa\ }
|\mathrm{both\ blank},W0,B1\mathrm{occupied}\rangle.
$$

All three states have on-site energy `2 Delta`. The hopping eigenvalues
are zero and `plus or minus sqrt(g^2+kappa^2)`. At duration

$$
t_{\mathrm{return}}=\frac{2\pi}{\sqrt{g^2+\kappa^2}}
$$

the defect returns exactly, up to a global phase. Starting instead from the
middle state, the occupied bath creates a discrepancy with probability
`g^2 t^2+O(t^4)` at short time, without another replica flip. Reverse
transport is part of the same physical law. Every finite version also has
unitary recurrences, so the guarantee has the declared finite horizon.

## 11. Fault work, entropy and general storage bounds

At the first fault after nominal acquisition, the ready triple and vacuum
drains have zero parity energy. An S flip excites both parities and has
mean work `2 Delta`. Either memory flip excites one and has mean work
Delta. Kinetic and writer-exchange expectations vanish in these initial
basis states. Later pulses can act on residual or returning excitations.
Their exact work is

$$
W_j=\operatorname{Tr}\rho\,(F_j^\dagger H F_j-H),
$$

with the on-site parity contribution
`Delta sum_i m_(j,i)(1-2<p_i>)` on ready states. A constant positive work per
later pulse is not assumed. The external pulse and timing resources remain
part of the experiment.

The two chains are initialized once and remain inside the system. Their
4608 pure spins, the blank memories, the weak writer coupling and its long
preparation time are declared resources. Global von Neumann entropy is
preserved by writing, correction and each specified unitary fault. It is
initially one bit for the fair source mixture. Random choices of faults
add classical uncertainty to an ensemble only when the corresponding
controller or choice record is excluded. Entropy transport into internal
degrees of freedom supplies no external entropy sink.

For comparison with larger alphabets, an exact three-replica single-fault
corrector for d classical labels has `3(d-1)` different single-substitution
inputs for each codeword, plus the undamaged input. If one common unitary
restores them to the same codeword from one environmental preparation,
unitarity requires their
environmental record states to be mutually orthogonal. The environment
therefore needs dimension at least `3d-2`, or `log2(3d-2)` bits of capacity.

| Alphabet d | Necessary exact syndrome-record dimension |
|---:|---:|
| 2 | 4 |
| 3 | 7 |
| 4 | 10 |
| 8 | 22 |
| 16 | 46 |
| 64 | 190 |

Two replicas have an ambiguity for an unrestricted unflagged substitution:
the same pair of distinct labels can arise from either historical label
with one changed site. This bound concerns that fault class. A flagged
erasure supplies additional information. The larger-alphabet statements
are necessary discrimination and storage bounds, and the numerical D4
probability certificate above applies to the specified binary model.

## 12. Verdict and clock scope

One fixed finite Hamiltonian acquires two classical records, recognizes the
declared single-site discrepancies through its supplied comparisons, and
repairs a sequence of up to four faults with the proven terminal historical
and resource scores. The same drains are reused and all backflow remains
in the model. The repair interaction begins upon a discrepancy through its
existing matrix elements.

The laboratory faults and terminal read are externally timed. The model
supplies neither their autonomous generation nor an internal physical
clock. Majority identification is valid in the certified separated-fault
class, and finite recurrences bound the lifetime. Recognition, correction,
local geometry and preparation are supplied standard-physics structures.
An emergent recognition rule, an additional GRUT law and a cosmological
selection principle remain separate, unearned verdicts.

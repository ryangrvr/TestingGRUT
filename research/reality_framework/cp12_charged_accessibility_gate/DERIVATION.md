# Derivation: charged cyclic access and collective probability

## A. Frozen physical space and target

There are N vertex fermions and M=binom(N,2) potential-edge fermions. Exactly Q graph edges are active. The total fermion number is M; every inactive edge has its edge fermion occupied and both half-leg fields zero. Consequently the active vertex/edge subsystem has Q fermions.

For an active edge e=(i,j), Gauss gives

\[
E_{ie}+E_{je}=1-n_{b_e}=h_e,\qquad
\sum_{e\ni i}E_{ie}=n_{a_i},
\qquad \sum_e h_e=n_v.
\]

The frozen Hamiltonian is

\[
H=\kappa Q+D+R+K,\quad
D=2m n_v+\frac{\epsilon_E}{2}\sum E^2,\quad
K=H_{\rm hop},
\]

with kappa=1/2, m=t=1, epsilon_E=1/4, Gamma=1/10. R is the original neutral-rewiring Laplacian with its full diagonal. No spatial weighting is introduced.

The normalized neutral ray is the equal coherent sum of all labelled Q-edge graphs with every edge fermion occupied, all vertex fermions absent, and E=0. It satisfies

\[
(D+R)\psi_0=0,\qquad \|K\psi_0\|^2=2Qt^2.
\]

Each neutral parent has 2Q distinct first-hop targets. Their graph, edge hole and occupied vertex identify a unique inverse parent, so no cross-parent cancellation occurs in this norm.

Define C(E) as the union of directed electric circuits: an edge is directed by its positive half-leg field and belongs to C only if its two fields have opposite nonzero signs and it participates in a directed cycle. Let e_C=|C(E)|.

For N>=5,Q=2N, the main candidate G_N requires connected four-regular g, e_C=Q, and n_v>=ceil(N/4). No particular host is selected by this definition. Its underlying cycle rank is Q-N+1=N+1. The upper-bound event F_N only requires e_C=Q and the same matter threshold; G_N is a subset.

## B. A bounded hopping operator despite unbounded fields

Fix a graph g and Fourier-transform the active electric rotors. At rotor angles theta, K is a number-conserving fermion quadratic form. Its one-particle matrix is

\[
t\begin{pmatrix}0&B(\theta)\\B(\theta)^\dagger&0\end{pmatrix},
\]

where B is N by Q and each column has two unit-modulus entries, with the frozen incidence signs. Therefore

\[
\operatorname{tr}(BB^\dagger)=2Q,\qquad
r_B\le r_*=\min(N,Q).
\]

If its nonzero singular values are s_j, the one-particle energies are +/- |t|s_j and zeros. On the Q-fermion active Fock sector the absolute energy is at most |t| sum_j s_j: filling the positive modes and enough zero modes attains that upper estimate in the unconstrained Fock space. Cauchy-Schwarz gives

\[
\boxed{\|K(\theta)\|\le |t|\sum_j s_j
\le |t|\sqrt{2Qr_*}=:b.}
\]

This gives ||K||<=b uniformly in g and theta. Restricting to the Gauss-law sector cannot increase the bound.

D is a nonnegative diagonal self-adjoint operator. R is bounded for fixed N,Q; a row-sum estimate gives ||R|| <= 2 Gamma Q(M-Q)/M. Thus the physical H is a bounded self-adjoint perturbation of D plus a constant, on the usual diagonal domain. The following argument does not require an electric cutoff or an analytic-vector assumption.

## C. Minimum electric distance and hopping order

Let

\[
\mathcal L(E)=\sum_{i,e}|E_{ie}|.
\]

Every through-flow occupied edge has half-leg L1 at least 2. Every through-flow hole has L1 at least 3, because its opposite-sign integer fields sum to 1. Every hole outside C has L1 at least 1. Hence, writing h_C for holes in C,

\[
\mathcal L\ge 2e_C+h_C+(n_v-h_C)
=\boxed{2e_C+n_v}.
\]

R preserves the entire electric tuple: it removes and inserts only zero-field edges. D is diagonal. Both preserve each L1 layer. One matter hop changes one integer field by +/-1, and therefore changes L1 by exactly one.

For a target P with e_C>=L and n_v>=k, set d=2L+k. It annihilates every L1 layer below d. Since H can raise L1 by at most one,

\[
PH^j\psi_0=0,\qquad 0\le j<d.
\]

This is a lower bound on the first possible order. Nonzero amplitudes at order d require an independent interference calculation.

Finite-support physical vectors belong to Dom(H^j) for every finite j: each application has finite support and finite fields. Strong finite-order Taylor expansion is therefore available even though convergence of an infinite Taylor series is not asserted.

## D. Collective probability bound

Set H0=kappa Q+D+R and K_I(s)=exp(iH0 s) K exp(-iH0 s). H0 preserves every L1 layer. Therefore K_I(s) is still nearest-layer, with norm at most b.

The interaction-picture propagator has a finite Duhamel expansion with unitary remainder:

\[
U_I(T,0)=\sum_{j=0}^{d-1}(-i)^j
\int_{0<s_1<\cdots<s_j<T}
K_I(s_j)\cdots K_I(s_1)\,ds
+\mathcal R_d,
\]

\[
\mathcal R_d=(-i)^d
\int_{0<s_1<\cdots<s_d<T}
U_I(T,s_d)K_I(s_d)\cdots K_I(s_1)\,ds.
\]

After multiplication by P exp(-iH0 T), every displayed term with j<d vanishes. The remaining leftmost propagator is unitary.

Because H0 psi0=kappa Q psi0,

\[
\|K_I(s_1)\psi_0\|=\|K\psi_0\|
=\sqrt{2Q}|t|=:\sigma.
\]

Thus, for d>=1,

\[
\|P e^{-iHT}\psi_0\|
\le \frac{\sigma b^{d-1}|T|^d}{d!},
\]

and

\[
\boxed{
p_P(T)\le\min\left\{1,\frac1{r_*}
\left[\frac{(b|T|)^d}{d!}\right]^2\right\}.
}
\]

For Q=2N and the full-circulation charged target, d_N=4N+ceil(N/4), b=2N|t|, r_*=N. Using d! >= (d/e)^d,

\[
p_{G_N}(T)\le
\min\left\{1,\frac1N
\left[\frac{2eN|tT|}{d_N}\right]^{2d_N}\right\}.
\]

For each fixed |tT|<4.25/(2e), the bracket is bounded below one at sufficiently large N, proving exponential suppression. Nothing follows at later times from this bound alone. The threshold is constant in N because the extensive hopping order and the extensive hopping norm compensate. This is not propagation through an intrinsic spatial metric.

## E. Exhaustive exact N=5 charged onset

At N=5,Q=10, g=K5 and R=0. The target threshold is n_v>=2. Any minimum-L1 target has n_v=2, L1=22, and every edge has the minimum through-flow magnitude: occupied edges carry (1,-1), holes carry (2,-1), up to endpoint reversal.

Replace each hole's 2 by 1 to obtain a tournament orientation. Let h_i be the number of holes with positive endpoint i. Its unit-field vertex sum is even (degree four), and Gauss says

\[
(\text{unit-field sum at }i)+h_i=n_{a_i}.
\]

With two holes, coincident positive endpoints would force every n_ai even, contradicting n_v=2. The positive endpoints are therefore distinct. Their vertices have n_ai=1, every other vertex has n_ai=0, and every unit-field sum vanishes. The tournament must be regular, with two positive and two negative half legs per vertex.

There are 24 regular tournaments on five labelled vertices. For each, the two holes can be any pair of its ten directed edges except the five pairs sharing a positive endpoint:

\[
\binom{10}{2}-5\binom22=40.
\]

A connected balanced directed graph has every edge on a directed cycle. Thus all 24 x 40=960 configurations qualify, and this exhausts the order-22 target support.

For each configuration, every shortest path is monotone in each electric half leg. The integer dynamic program records whether each leg has been visited zero, once, or (on a hole's positive leg) twice, together with fermion bits. At each allowed hop it multiplies by the frozen endpoint sign and the fermionic creation/annihilation sign, and adds amplitudes. No probabilities or absolute path counts are substituted.

The eight permutation orbits each have 120 configurations. The amplitudes have absolute values

\[
a=49583226880,\qquad b_0=40677670912,
\]

each in four orbits. The full coherent squared norm is

\[
480(a^2+b_0^2)=1974321263356092214149120.
\]

All 22 events are hops, so neither the diagonal terms nor rewiring can contribute to this first nonzero amplitude. Dividing by (22!)^2 gives

\[
p_{G_5}(T)=
\frac{1326415312}{848779513441154588027109375}T^{44}
+O(T^{46}).
\]

Reality of H, psi0 and P_G makes the probability even in T, eliminating the possible order-45 cross term. This statement uses finite-order differentiability, not a claimed globally convergent Taylor series.

The reverse dynamic program removes electric increments from the target. It uses a separate fermion-sign rule based on occupied modes strictly between the source and destination. Every orbit agrees with the forward program. All 960 labelled forward amplitudes were additionally calculated directly; their squared sum independently confirms the orbit aggregation. One identical labelled target was checked in Python and JavaScript.

## F. Actual cyclic support, rather than energy-moment feasibility

Let Pi_cyc be the orthogonal projector onto the closed cyclic subspace generated by the spectral functions of H acting on psi0. Since psi0 is in Dom(H^22), (H-5)^22 psi0 belongs to that subspace. The order obstruction gives

\[
P_G(H-5)^{22}\psi_0=P_GH^{22}\psi_0\ne0.
\]

Thus P_G Pi_cyc is nonzero for the **actual original ray**. This goes beyond constructing a different state with the same first two energy moments.

For an explicit conservative lower bound, vectors reached in j applications have L1<=j, n_v<=min(5,j), and sum E^2<=j^2. At N=5,Q=10, ||K||<=10 and R=0, so

\[
\|(H-5)\Pi_{\mathcal L\le j}\|
\le 10+2\min(5,j)+j^2/8=:c_j.
\]

The first step has exact norm sqrt(20). Therefore

\[
\|(H-5)^{22}\psi_0\|^2
\le 20\prod_{j=1}^{21}c_j^2=:C^2,
\]

and

\[
\|P_G\Pi_{\rm cyc}\|^2
\ge \frac{1974321263356092214149120}{C^2}
=3.272338151350656\ldots\times10^{-42}.
\]

The rational value is recorded in EXACT_ONSET.json. This is an operator-overlap bound; it does not bound p_G(T) from below uniformly in time and does not establish appreciable spectral weight or a thermodynamic phase.

## G. Full-H residual control of an initial-weighted spectral envelope

Take any finite collection of physical trial vectors v_a in Dom(H), real trial frequencies lambda_a and coefficients c_a. Define

\[
z(T)=\sum_a c_a e^{-i\lambda_aT}v_a,\quad
e_0=\|\psi_0-z(0)\|,\quad
r_a=\|(H-\lambda_a)v_a\|.
\]

Variation of constants gives the exact estimate

\[
\|e^{-iHT}\psi_0-z(T)\|
\le e_0+|T|\sum_a|c_a|r_a=:\epsilon(T).
\]

No spectral completeness, degeneracy assumption, or infinite-dimensional diagonalization is needed. For a target projector,

\[
\|P_Gz(T)\|^2
=\sum_{a,b}e^{i(\lambda_a-\lambda_b)T}
c_a^*c_b\langle v_a|P_G|v_b\rangle
\le A,
\]

where A is the entrywise absolute sum of the **initial-weighted target Gram matrix**. Consequently

\[
p_G(T)\le\min\{1,(\sqrt A+\epsilon(T))^2\},
\]

and for a computed trial probability q(T),

\[
(\sqrt q-\epsilon)_+^2
\le p_G(T)\le\min\{1,(\sqrt q+\epsilon)^2\}.
\]

If P retains a finite physical space and B=(1-P)HP contains all outgoing transitions, then

\[
r_a^2=
\|(PHP-\lambda_a)v_a\|^2+\|Bv_a\|^2.
\]

This uses the full operator. An all-time bound on the projected trial signal must not be called an all-time bound on the full physical evolution after dropping epsilon(T).

For the N=4,Q=4 signed S4 benchmark, all 577 numerical eigenvectors of PHP are used and c_a=<v_a,psi0>. The actual ray enters every weight. Every full-law outgoing transition is included in B. The values and numerical qualifications are in NUMERICAL_AUDIT.json.

## H. Finite-time capacity check and limitations

The N=5,Q=10 budget32 retained model contains the entire minimum-order target support (sum E^2=26), but not every possible charged target: n_v=5 full circulation already needs sum E^2>=35. The parent leakage method bounds the norm difference between retained evolution and the full-H state. The projector-contraction interval in Section G then includes omitted target support.

At T=1 the executed retained probability is 2.1256306851140608e-19 and the full-state error upper estimate is 8.539230599708589e-7, giving the full upper estimate 7.299721998951689e-13. Its lower bound is zero. A positive retained probability is not promoted to a resolved positive full-law numerical response.

The exact onset and Duhamel bounds are analytical, with exact integer checks of onset. The spectral and leakage estimates are numerical evaluations of analytical inequalities; no formal machine-rounding enclosure has been constructed. The packaged Python rerun is pending after the container disconnected. There is no claim of external proof certification, autonomous preparation, extended persistence, dimension selection, or a spatial causal cone.

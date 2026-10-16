# S6: complete cross-history receiver calculus

This note states exact-arithmetic identities for the unchanged S4/S5
formation Hamiltonian. Numerical validation, computational efficiency,
and all-size formation are separate obligations.

## Orthogonal histories, distinct receiver frames

Write H=A+B as in S5. For each type word u of length k choose an
orthonormal L_N by r_u frame W_u. Let F_k(W_u) be its normalized bosonic
Fock isometry. The numerical isometry and physical projector are

    V = direct_sum_u |u><u| tensor F_|u|(W_u),  P=V V^dagger.

History labels stay orthogonal. Different receiver frames are generally
neither equal nor mutually orthogonal. Every coefficient vector retains
the complete ordered history and all occupations with exactly |u| bosons.
The vacuum history has a one-dimensional coefficient space.

For any rectangular contraction G, Gamma_k(G) is its symmetric tensor
power in normalized occupation bases. In particular

    F_k(W_v)^dagger F_k(W_u)=Gamma_k(W_v^dagger W_u).

This includes every permanent, multiplicity and phase. It is not the
identity when the two frames differ. For child v=ua define
lambda_v=(g/N)sqrt(n_a-k_a(u)) and c_v=W_v^dagger c, where c is physical
contact. Then the complete projected append block is

    (V^dagger B V)_(v,u)
      = -lambda_v a^dagger(c_v) Gamma_k(W_v^dagger W_u).

The removal block is exactly its adjoint. Equivalently it acts as
Gamma_k(W_u^dagger W_v) a(c_v). No order history, reverse or field
overlap has disappeared. PBP remains Hermitian and has norm <=g.

## Full residual, including sibling interference

For a coefficient vector x, let phi_u=F_|u|(W_u)x_u. Its full physical
emission action in target history u is

    y_u= -lambda_u a^dagger(c) phi_parent(u)
         -sum_(v child of u) lambda_v a(c) phi_v.

The first term is absent at the root. Projection gives the preceding
append/pop blocks. Since histories are orthogonal,

    ||R x||^2 = sum_u ||y_u||^2 - ||V^dagger B V x||^2,
    R=(1-P) B V.

All terms in ||y_u||^2 are computable without physical Fock enumeration.
For a parent vector p, creation contributes

    ||a^dagger(c) phi_p||^2 = ||x_p||^2+||a(c_p)x_p||^2.

For children v,w put z_v=a(c_v)x_v. Their cross term is

    <a(c)phi_v,a(c)phi_w>
       = <z_v,Gamma_|u|(W_v^dagger W_w) z_w>.

Creation/removal interference is

    <a^dagger(c)phi_p,a(c)phi_v>
       = <x_p,Gamma_(|u|-1)(W_p^dagger W_v) a(c_v)^2 x_v>.

Multiply by the actual lambda coefficients and include conjugates.
Discarding these cross terms would change the residual, even though
the original histories themselves are orthogonal. Their Hamiltonian
images can occupy the same target history. Near cancellation requires
careful arithmetic; these identities alone are not interval bounds.

## Projected-emission error with all future returns

Hold frames fixed during a B step h>=0. Put C=V^dagger B V, ||C||<=g,
and start with normalized numerical x. The exact Duhamel identity gives

    ||exp(-iBh)Vx - V exp(-iCh)x||
       <= integral_0^h ||R exp(-iCs)x|| ds.

For any integer q>=1, expansion of exp(-iCs) through q-1 yields

    rho_h(x) = sum_(l=0)^(q-1) h^(l+1)/(l+1)! ||R C^l x||
               + exp(gh)(gh)^(q+1)/(q+1)!.

The integral is at most rho_h. Indeed ||R||<=g and the exponential
Taylor remainder at s is <=exp(gs)(gs)^q/q!; integrating pays the
displayed tail. S6 uses q=3. This residual is state dependent but the
uncomputed tail is fully bounded; no residual quadrature is guessed.

A degree-p polynomial for exp(-iCh)x adds at most
exp(gh)(gh)^(p+1)/(p+1)!. Normalize the resulting global vector once
and charge its observed norm change. S6 uses p=8.

Between B steps, A is exact: each word gets its own energy phase and
each W_u gets the same physical one-particle free unitary. All L_N
frequencies remain represented. S5's physical-H commutator constant
C_Strang(g) is unchanged. For n steps at T=n h,

    full-state error <= C_Strang T^3/n^2 + sum_steps rho_h
                        + sum orbital-projection changes
                        + all exponential/normalization/frame changes.

Physical unitarity telescopes these errors even when frames are chosen
adaptively from prior numerical states. Every omitted amplitude's entire
future action, including coherent returns, is paid. This is a bound on
the original physical model, not locality obtained by removing paths.

## History-conditioned mode projection

An orbital-vacuum projection in a particular history discards global
probability q_u measured with that history's unnormalized amplitude.
For simultaneous projections in distinct histories, q=sum_u q_u.
Global projection followed by one global normalization changes state by

    sqrt(2q/[1+sqrt(1-q)]).

Word amplitudes are never separately normalized. Natural-orbital
selection is only a numerical choice of subspace; it cannot select a
favorable physical root, field, geometry or observation time.

## Cost, without a polynomial-complexity claim

With D(r,k)=binomial(r+k-1,k), coefficient storage is

    sum_u D(r_u,|u|),  frame storage L_N sum_(u nonempty) r_u.

Direct cross-Fock matrices cost D(r_v,k)D(r_u,k) complex entries per
needed pair. Besides parent/child edges, full residuals need sibling
pairs and grandparent/grandchild pairs. Factorial history growth is
still present for all distinct detunings.

To avoid large dense blocks, factor G=U S V^dagger. Apply the exact
Fock lifts of the two unitaries using two-mode rotations and multiply
each occupation by the product of its singular values. Rectangular
zero singular directions are projected exactly as part of Gamma(G).
This is algebraic evaluation of the cross-frame contraction, not an
additional physical truncation. Rotation groups cost at most polynomial
in rank times D(r,k); they do not remove D's k-dependence or history
growth. Neither uniformly small ranks nor scalable formation follows.

An actual solver can fail its fixed error/resource budget while these
transition/error identities remain valid. Such failure must be recorded
separately from failure of the physical formation mechanism.

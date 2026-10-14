# Observable error control for the CP-12 source response

This packet tightens the numerical error control of the benchmark at `4570ed1c2ac2f8dab3f49ddb99c6224865eb1527`. It does not change H, the electric-square budget, the neutral preparation, the waiting time, or the source intervention. The argument is a general projection-error identity, specialized to the declared CP-12 observables. Numerical bounds are evaluated in ordinary floating point, not interval arithmetic.

## Forward residuals

Let P project onto the retained physical space and let H_P=PHP and B=(1−P)HP. Work within the exact signed source-fixing symmetry quotient. It is a physical subspace of the full theory; the projection here is the electric-budget/reachability projection inside that sector. Every full-law outgoing transition is present in B. All frozen rewiring diagonals are retained in H_P.

For any finite retained vector a, Duhamel bounds its forward error by

    ε(a,t) = ∫₀ᵗ ||B exp(−iH_P s)a|| ds.

The parent midpoint Taylor upper sum evaluates this integral from above analytically before floating-point evaluation. A reversed segment uses −H_P; its leakage matrix is −B, with identical norm. No thermal or bounded-full-Hamiltonian approximation is introduced.

Now consider two paths from a normalized retained initial vector ψ:

    x_exact = U(τ)U(w)ψ,
    y_exact = U(τ)S U(w)ψ,

where S is a contraction commuting with P. For the response, S=n₀ is the vertex source occupation, not a unitary. It arises as the derivative of the unitary intervention exp(−iηn₀). Denote their projected counterparts by x and y.

If ε_w is the pre-source propagation bound and ε_x,lag and ε_y,lag are the two post-source bounds, then

    ||x_exact−x|| ≤ e_x = ε_w+ε_x,lag,
    ||y_exact−y|| ≤ e_y = ε_w+ε_y,lag.

The insertion S cannot increase the pre-source error. The previous benchmark bounded the response directly by a quantity linear in e_x+e_y. That is safe but very loose for these observables.

## Use the detector's compatibility with P

Let O be Hermitian, with [O,P]=0 and a known GLOBAL norm upper bound o. For the detector averages used here, o=1 on the entire unregulated physical Hilbert space. Using only the largest value of O over retained basis states would not in general establish that global bound.

Put δ_x=x_exact−x and δ_y=y_exact−y. Exactly,

    ⟨x_exact,O y_exact⟩−⟨x,O y⟩
      = ⟨x,O δ_y⟩ + ⟨δ_x,O y⟩ + ⟨δ_x,O δ_y⟩.

The third term is bounded by o e_x e_y. For the first two terms, use a backward residual rather than merely the norm of x or y.

Take the first term and set z=O x. Since [O,P]=0, z is retained. The forward error δ_y is an integral of omitted-space residuals B v(s), followed by the remaining exact unitary segments and, when appropriate, S. Move those remaining operators onto z. Their projected backward evolution remains in P, hence has zero inner product with B v(s). Only the ERROR of that backward evolution contributes.

Let d_y bound the entire backward-path error from O x, first through a reversed lag segment, then S†, then a reversed waiting segment. Its cumulative bound dominates the backward error to every intermediate time because every insertion is a contraction and the residual upper integrals are nonnegative. Therefore

    |⟨x,O δ_y⟩| ≤ e_y d_y.

Similarly, let d_x bound backward evolution from O y over the uninserted total duration w+τ. Then

    |⟨δ_x,O y⟩| ≤ e_x d_x.

Combining the three terms proves

    |⟨x_exact,O y_exact⟩−⟨x,O y⟩|
       ≤ e_y d_y + e_x d_x + o e_x e_y.

There is no claim that this always behaves quadratically in a single small parameter: detector-dependent backward errors can be large. Its advantage must be computed for each detector. It is valid even when O does not commute with H or with S; implementation here uses physical diagonal source/detector operators.

## Response bound

The source-kick linear susceptibility is

    χ(O) = 2 Im ⟨x_exact,O y_exact⟩.

Thus the projected value χ_P=2 Im ⟨x,O y⟩ has absolute error at most

    E_O = 2(e_y d_y + e_x d_x + e_x e_y),

using the global detector norm upper bound one. In the code,

    d_y ≤ ε_backward(O x, lag)
           + ε_backward(S† exp(+iH_P lag) O x, wait),
    d_x ≤ ε_backward(O y, wait+lag).

The same parent quadrature upper-sum calculation is used in all four directions. This is a finite-time error argument. No fixed-N Taylor coefficient is substituted into a large-time evolution.

## Why the assumption matters

If O couples retained and omitted spaces, O x is not retained and the crucial residual orthogonality fails. A two-state counterexample makes this explicit: H=σ_x, P=|0⟩⟨0|, ψ=|0⟩ and O=σ_y. At t=0.1 the actual expectation differs from the projected one by |sin 0.2|≈0.1987. Incorrectly using the block-diagonal argument with P O P=0 would give only t²=0.01. The checker records this as a required applicability rejection, not as a valid bound.

## Detector interpretation and a metric-independent control

The first two detector averages are unchanged from the parent:

    O_r = (1/(N−1)) Σ_{j≠0} 1[d_G(0,j)=r] n_j,    r=1,2.

They are gauge invariant and commute with the electric projector. They use the FINAL graph distance, so a response can contain graph redistribution. Resolving their signs does not isolate transport at a previously fixed distance or establish a causal cone.

At Q=2N and N=5,6, every allowed graph is connected and has diameter at most two. Consequently the following operator identity is exact on the entire physical space:

    O₁+O₂ = O_other = (1/(N−1)) Σ_{j≠0} n_j.

The code also evaluates O_other directly, with its own backward error control, rather than inferring its uncertainty by subtracting bins. O_other has no graph-distance projector. A nonzero response in this channel cannot be explained solely by reassigning the same density among final-distance bins.

Matter hopping is necessary for every such response. With t_hop=0, the remaining mass/electric diagonal and rewiring operator commute with n₀, and both O_r and O_other commute with n₀. Therefore [U† O U,n₀]=0 and the source-kick susceptibility vanishes for every state. A nonzero O_other response witnesses a coupled matter-density effect. This is an operator attribution control, not a replacement physical theory or evidence that matter has selected spatial topology.

The source is local in the model's vertex occupation algebra, and the intervention is gauge invariant. Neither locality of that algebra nor these tiny, diameter-two graphs supplies a demonstrated extended spatial locality or a relativistic signal cone. The topology-neutral initial pure ray remains a supplied boundary condition.

# B2: detector reconstruction from cell averages

This is one candidate correction to the **assumed B1 detector dictionary**. It does not derive that dictionary from the quantum Hamiltonian. No history probability, target discrepancy, or fitted parameter enters the weights.

Let a=L/(N+1), I_j=[(j+1)a,(j+2)a], j=0,...,N-2, and y_j be a cell average. The two boundary values are y(0)=y(L)=1/4. The normalized Gaussian f, interval L=5, and sigma=0.8 are inherited. Reconstruct a piecewise quadratic R_a y:

* On each interior cell j=1,...,N-3, set xi=(z-(j+3/2)a)/a. The coefficients c of q(xi)=c0+c1 xi+c2 xi^2 solve M c=(y_{j-1},y_j,y_{j+1}), with rows M_r=(1,r,r^2+1/12), r=-1,0,1.
* On [0,2a], set xi=z/a and solve M_edge c=(y(0),y_0,y_1), with rows (1,0,0), (1,3/2,7/3), (1,5/2,19/3). Reflect this construction onto [L-2a,L].

The intervals cover [0,L] once, up to endpoints. Each reconstructed cell has its specified average. The first/last polynomial also satisfies the supplied boundary value. Adjacent reconstructed polynomials need not be continuous. Integrating f times these polynomials defines coefficients w_j,b_left,b_right:

    F_B2 = sum_j w_j y_j + b_left y(0) + b_right y(L).

For each stencil, compute moments m_p=integral f(z) xi^p dz and solve M^T coefficients=m. This fixes all weights by geometry and moments. Constant, linear, and quadratic profiles are reproduced exactly. In particular sum(w)+b_left+b_right=1, and reflection symmetry follows from symmetric f and mirrored stencils.

## Classical consistency and its limitation

For a scalar C^3 profile with the supplied boundary values, Taylor's theorem and bounded inverse stencil matrices give an O(a^3) uniform reconstruction error. An explicit conservative bound is C a^3 ||y'''||_infinity, where

    C_int = ( (1/2)^3 + Lambda_int (3/2)^3 ) / 6,
    C_edge = ( 2^3 + Lambda_edge 3^3 ) / 6,
    C = max(C_int,C_edge),

and Lambda is the supremum of the sum of the absolute entries of (1,xi,xi^2) M^{-1}, on xi in [-1/2,1/2] or [0,2], respectively. These are finite constants independent of a. Averaging the Taylor remainder does not enlarge its supremum; interpolation is exact on its quadratic Taylor part. Since f is nonnegative with integral one, the same bound controls the integrated scalar functional error.

A quantum field is an operator-valued distribution. This smooth-profile argument is **not** a quantum operator error bound or a proof of convergence of instruments. The actual free-vacuum Gaussian histories are evaluated separately, retaining their ordered complex covariances.

## Quantum implementation and matched reference

Keep the inherited Hamiltonian, neutral sector, ghost-endpoint geometry, and mass shift m_lat=m_R-a/8 unchanged. Gauss law gives L_j=S_j+1/4+(-1)^j/4, where S_j=(1/2)sum_{i<=j} Z_i. Under the B1 cell-average interpretation y_j=L_j-(-1)^j/4:

    F_lat,B2 = sum_j w_j (L_j-(-1)^j/4) + (b_left+b_right)/4
             = 1/4 + sum_j w_j S_j.

The matched continuum functional is integral f_B2(z) Phi(z) dz +(b_left+b_right)/4, with f_B2=w_j/a on I_j and zero in the two boundary gaps. This exactly represents the reconstructed functional of cell averages. The target remains integral f(z) Phi(z) dz. Neither target nor Hamiltonian is retuned.

For sine modes k_n=n pi/L, the matched coefficient is computed analytically:

    f_n = sqrt(2/L) sum_j (w_j/a) [cos(k_n lo_j)-cos(k_n hi_j)] / k_n.

Use the inherited ordered Gaussian-word formula with omega_n=sqrt(k_n^2+1/pi) and h_n^2=f_n^2/(2 pi omega_n). The total variation B of the zero-extended step smearing bounds the omitted variance by B^2 L^2/(2 pi^4 M^2). The inherited conservative r-event history bound is 2^r r^2 coupling^2 times this variance bound. These analytic bounds are evaluated in floating point, alongside mode doubling and quadrature diagnostics.

## What the error budget means

For every registered history, compare lattice to matched reference and matched reference to target separately. Their sum plus the frozen numerical reserve must be at most 0.001. A small direct lattice-to-target discrepancy cannot compensate for either component exceeding its assigned budget. This prevents the B1 cancellation from being counted as success. Solver residuals and refinement are numerical evidence, not interval-certified error estimates for the full physics.

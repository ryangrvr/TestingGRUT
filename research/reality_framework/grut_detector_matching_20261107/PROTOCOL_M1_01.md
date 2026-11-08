# M1-01 — independent massive interval reference and bounded history comparison

Freeze this document, card, code and analysis before inspecting any M1-01 histories. M0-01 and M0-02 are known misses and remain so. The latter motivates using its analytic ghost-endpoint improvement here; this is not a blind selection of that dictionary. The new calculation is a reference-construction and finite-resolution test, not a claim that the massless bridge has passed.

## Physical experiment

Use units hbar=c=g=1, interval L=5, theta=0, neutral sector and Dirichlet values phi(0)=phi(L)=sqrt(pi)/4. Fix renormalized fermion mass m/g=0.5. Prepare the interacting ground state independently in each formulation. The target observable is F=integral f(z) phi(z)/sqrt(pi) dz, with the same normalized Gaussian of width 0.8 as M0-02. The instruments are K_s=(I+s exp(i c F))/2, c=2, s=+/-1. Retain measurement disturbance. Evaluate all outcomes at times (0), (0,1), (0,2), (0,1,2). No new coupling, phase, mean, spectrum or state is fitted.

The lattice uses the reviewed Gauss-law Hamiltonian, incoming flux zero, even neutral chain, a=5/(N+1), m_lat=0.5-a/8, sites (j+1)a and link midpoints (j+1.5)a. The primary observable is the normalized Gaussian-weighted raw link field. The finite-interval mass prescription and endpoint map remain assumptions with residual cutoff error; M0-02 did not certify them to 0.1%.

## Independent continuum Hamiltonian

Adopt the standard bosonized Schwinger Hamiltonian, with infinite-line normal ordering at mu=g/sqrt(pi):

$$
H=H_0-\frac{m g e^\gamma}{2\pi^{3/2}}\int_0^L dz\,:\cos(2\sqrt\pi\phi(z)):_\infty,
\qquad H_0=\sum_{n\ge1}\Omega_n a_n^\dagger a_n+\text{constant}.
$$

Here Omega_n=sqrt((n*pi/L)^2+1/pi), u_n=sqrt(2/L) sin(n*pi*z/L), and

$$
\phi(z)=\bar\phi(z)+\sum_n\frac{u_n(z)}{\sqrt{2\Omega_n}}(a_n+a_n^\dagger),\qquad
\bar\phi(z)=\frac{\sqrt\pi\cosh[\mu(z-L/2)]}{4\cosh(\mu L/2)}.
$$

The shift diagonalizes H0; it is not imposed as the interacting expectation value. The cosine changes that expectation through its ground state. The common free vacuum constant is removed. Absolute interacting vacuum energy is not a renormalized observable in this calculation; only gaps and normalized histories are compared. No uncalculated Hamiltonian-truncation counterterms are silently included.

The starting action, normal-order convention and interval boundary dictionary are from [Okuda, equations (11), (13), (21)–(23), (66)](https://arxiv.org/abs/2210.00297). The finite-cutoff implementation below is independently constructed, not that paper's DMRG calculation. The bulk mass-improvement prescription is from [Dempsey et al.](https://arxiv.org/abs/2206.05308); its open-boundary limitations remain explicit.

## Normal ordering on an interval

Let beta=2 sqrt(pi), psi=phi-bar(phi), and normal order psi in the Dirichlet vacuum. Contraction subtraction gives

$$
:e^{i\beta\phi(z)}:_\infty=e^{i\beta\bar\phi(z)}\lambda(z):e^{i\beta\psi(z)}:_D,
\quad \log\lambda=2\pi[C_\infty(0)-C_D(z,z)].
$$

Using C_infty(r)=K0(mu |r|)/(2*pi) and Dirichlet images yields

$$
\log\lambda(z)=\sum_{j\in\mathbb Z} K_0(2\mu|z+jL|)-2\sum_{j\ge1}K_0(2\mu jL).
$$

This is checked against the independent integral representation in Okuda's equation (22). Near a boundary lambda grows as exp(-gamma)/(mu z); the phase beta*bar(phi) approaches pi/2, so the vacuum cosine vanishes linearly. Nontrivial normal-ordered matrix elements also have their Dirichlet sine factors. The finite matrix integrands have removable endpoint limits. Gauss–Legendre nodes remain strictly inside the interval; doubling quadrature checks resolution, not ultraviolet Hamiltonian convergence.

For a mode with real coefficient ell,

$$
\langle n|e^{i\ell a^\dagger}e^{i\ell a}|m\rangle=
\sqrt{\frac{\min(n,m)!}{\max(n,m)!}}
(i\ell)^{|n-m|}L_{\min(n,m)}^{|n-m|}(\ell^2).
$$

Multiply these factors across modes and take the real part after multiplying by exp(i beta bar(phi)). This builds the cosine without exponentiating a truncated field operator. Check its matrix elements against an independently exponentiated 80-level oscillator. Check image count 8 versus 16 and the integral formula at z=0.02,0.2,1,2.5,4.8.

## State and instrument cutoff

P_E retains every free occupation tuple with sum Omega_n n_n <= E. Use E=4,6,8,10,12; dimensions counted before production are 21,72,207,534,1280. The interacting Hamiltonian is P_E H P_E, diagonalized as a real symmetric matrix. Its lowest eigenvector is the preparation. This approximation requires its own refinement; agreement at m=0 does not bound interacting errors.

Compress the actual Weyl instrument, W_E=P_E exp(i c F) P_E. With h_n=integral f u_n / sqrt(2*pi*Omega_n), its matrix is

$$
(W_E)_{\boldsymbol n\boldsymbol m}=
e^{ic\bar F-c^2\sum_{k\ge1}h_k^2/2}
\prod_k\langle n_k|e^{ic h_k a_k^\dagger}e^{ic h_k a_k}|m_k\rangle.
$$

Higher unoccupied modes contribute their vacuum factor. The infinite sum uses the same independently computed 8192-mode reference and its analytic tail bound. W_E is a contraction, not generally unitary. Accordingly

$$
\sum_{s=\pm}K_{s,E}^\dagger K_{s,E}=\frac{I+W_E^\dagger W_E}{2},\qquad
K_{o,E}^\dagger K_{o,E}=\frac{I-W_E^\dagger W_E}{2}.
$$

Retained branches evolve with P_E H P_E. The overflow outcome enters an absorbing classical register; it is never renormalized away or turned into a physical detector outcome. Compare distributions on all binary strings plus an overflow symbol, with zero reference/lattice overflow. Check positivity of the defect matrix, total probability, and cutoff stability. Overflow quantifies this instrument compression only: it does not bound omitted Hamiltonian dynamics, coherent returns or a complete continuum error. No simultaneous-cutoff theorem is asserted.

## Frozen tests and resources

The executable thresholds are in [protocol_m1_01.json](protocol_m1_01.json).

- Verification: vertex error <1e-11; normal-factor disagreement <1e-9; Hermiticity <1e-10; eigen residual <1e-9; defect minimum >=-1e-10; probability normalization error <1e-9.
- Free verification at E=6,10: E=10 histories, including overflow, must be within TV 1e-6 of the independent Gaussian continuum reference.
- Massive reference: finest E=10-to-12 change <=1e-4 for every history; quadrature 128-to-256 at E=10 <=1e-6; maximum finest overflow <=1e-5. These are necessary numerical gates, not certified absolute error bounds.
- Lattice: N=10,14,18,20; tighter N=14 with split propagation must change histories <1e-8. Compare N=20 against E=12 with TV <=0.001 and N=18-to-20 change <=0.00025. Report all protocols, even if a prerequisite fails.
- One BLAS thread, maximum dimension 1300, 300 seconds per point, 1800 seconds for all subprocesses. No higher cutoff, mass fit, new endpoint shift or extrapolated substitute in this registration. Preserve failures and caps.

The comparison remains unresolved if any requisite numerical gate fails. A persistent mismatch with individually resolved approximations would motivate a dictionary/implementation investigation before a physical claim. Even a numerical pass would establish only these ideal finite histories in the adopted model and boundary sector. Extended instantaneous detectors have no causal implementation error bound here. Neither new GRUT dynamics nor empirical identification follows.

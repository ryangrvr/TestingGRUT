# Photon bridge clarification and selection verdict

7 November 2026 · Amendments to [PHOTON_MATCHING.md](PHOTON_MATCHING.md)

This derivation answers the occupation-cutoff, uniformity and locality questions for the existing free photon construction. It specifies a missing cutoff convention, proves exact axial matching with the completed instruments, and adds a restricted commutator result. The current GRUT postulates do **not** select plaquette-only dynamics. Both members of the comparison remain admitted models.

The analytic arguments and verification program are by the same author as the bridge. This is a response to the review, with reproducible checks, rather than independent reproduction.

## 1. Fix the regulator before comparing instruments

The normal-mode occupation projector is part of the finite model. It is not determined by the Hamiltonian eigenvalues alone. In a degenerate oscillator sector, a product cutoff in one mode basis generally differs from a product cutoff in a rotated basis.

Use the following convention for both Hamiltonians.

1. Fix the same odd spatial lattice, midpoint Fourier convention, representative of each nonzero pair of opposite momenta, and real cosine/sine quadratures.
2. Never mix different momentum pairs when diagonalizing the quadratic Hamiltonian, including across accidental frequency degeneracies.
3. For an axial wavevector, fix the same two coordinate-transverse polarization vectors in both models. Resolve their signs and ordering by the same coordinate convention. Do not rotate this axial basis independently between models.
4. Diagonalize the off-axis transverse blocks separately. Any residual choices there must be recorded but do not affect the axial theorem.
5. Cut each resulting real oscillator at the same occupations \(0,\ldots,K\).

This fixes the finite models needed for the comparison. It does not guarantee that the auxiliary product cutoff respects every spatial lattice symmetry. In particular, a symmetry mixing degenerate polarizations can fail to preserve a rectangular occupation region. Symmetry claims for the uncut spatial Hamiltonian must not be silently assigned to its finite occupation compression.

Let \(X\) denote the common collection of axial real oscillators and \(R\) the remaining ones. The untruncated Hilbert spaces, Hamiltonians and cutoffs factor as

\[
\begin{aligned}
\mathcal H_\eta&=\mathcal H_X\otimes\mathcal H_{R,\eta},&
H_\eta&=H_X\otimes1+1\otimes H_{R,\eta},\\
\Omega_\eta&=\Omega_X\otimes\Omega_{R,\eta},&
P_{K,\eta}&=P_{K,X}\otimes P_{K,R,\eta}.
\end{aligned}
\tag{R1}
\]

The factors in \(X\), including their canonical normalization and vacuum, agree because the added magnetic term vanishes on every axial mode. The spectator Hamiltonians and spectator vacua need not agree under the physical canonical-field identification.

An axial-only probe is \(W_\eta=W_X\otimes1\), so its compressed ordinary Kraus operators, on the finite range, are

\[
\widehat K_{s,\eta}=k_{s,X}\otimes1_{R,\eta},\qquad
k_{s,X}=P_{K,X}\frac{1+s e^{i\varphi}W_X}{2}P_{K,X}.
\]

Writing \(D_X=1_X-\sum_{s=\pm}k_{s,X}^\dagger k_{s,X}\), the completed overflow operator is exactly

\[
\widehat K_{{\rm ov},\eta}
  =\sqrt{D_X\otimes1_{R,\eta}}
  =\sqrt{D_X}\otimes1_{R,\eta}.
\tag{R2}
\]

Free evolution also factorizes. Every history branch, including branches after one or several overflows, therefore factors into the same axial branch and a spectator unitary. Taking its squared norm proves

\[
p^{(a,K)}_{\pi,\eta}(s_1,\ldots,s_m)
 =p^{(a,K)}_{\pi,0}(s_1,\ldots,s_m),
\qquad s_j\in\{+,-,{\rm ov}\}.
\tag{R3}
\]

This equality holds at every cutoff, for axial-only protocols and the stated vacua. More generally it needs only identical initial axial reduced states. It does not assert equality of the full finite instruments for off-axis probes, arbitrary new readouts, or independent choices of axial cutoff basis. Superpositions of different axial wavevectors are allowed. Such probes need not have localized spatial support.

### A control that fails without the convention

Take two degenerate oscillators with vacuum and \(K=1\). Displace only the first mode by a coherent amplitude \(\alpha\), with \(\mu=|\alpha|^2\). The probability that the displacement leaves the rectangular cutoff is

\[
\ell_{\rm original}=1-e^{-\mu}(1+\mu).
\]

In the equally valid diagonal basis \(b_\pm=(b_1\pm b_2)/\sqrt2\), the two amplitudes each have squared modulus \(\mu/2\). Cutting each rotated oscillator at \(K=1\) gives

\[
\ell_{\rm rotated}=1-e^{-\mu}(1+\mu/2)^2.
\tag{R4}
\]

For the completed two-outcome instrument of (P2), the first overflow probability is \(\ell/2\). At \(\mu=1\), the two probabilities are approximately \(0.13212056\) and \(0.08613563\). The Hamiltonian, vacuum and uncompressed probe were unchanged. This disproves basis-independent finite matching and explains why (R1) must be specified.

Compression means taking matrix elements of the **uncut** Weyl operator. Exponentiating a compressed field, \(e^{iPFP}\), defines a different instrument and generally removes the intended overflow. The replay uses exact oscillator displacement matrix elements, not that substitution.

## 2. Constants and quantifiers in P4

The claim concerns the prepared vacuum and a fixed finite list of protocols. For a protocol \(\pi\) fix its length \(m\), clock readings \(0\le t_j\le T\), real smooth torus smearings \(f_j=(f_{E,j},f_{B,j})\), probe strengths included in those smearings, and phases \(\varphi_j\). Fix \(0<c<\infty\), \(0<L<\infty\), and a finite \(\eta\ge0\). Fix all midpoint and transverse-plane identification conventions. No smearing is narrowed while taking the limits.

### Occupation error at fixed spacing

In the Heisenberg picture write \(W_j=D(z_j)\), where \(z_j=z_j(a,\eta,L,c,f_j,t_j)\) is the displacement vector of all independent real oscillators, in the convention of Section 1. Its norm obeys

\[
\|z_j\|^2=\langle\Omega_\eta|F_j(t_j)^2|\Omega_\eta\rangle.
\]

It is independent of the free time \(t_j\). Put

\[
d_j=\sum_{\ell=1}^j\|z_\ell\|,\qquad
B_\pi=\sum_{j=1}^m d_j,\qquad
C_{\rm occ}(\pi,a,\eta,L,c)=2^{m+1}B_\pi.
\tag{R5}
\]

Each uncut branch after \(j\) instruments is a sum of coherent states with coefficient absolute values summing to at most one. Their displacement norms are bounded by \(d_j\). Since a state outside the product cutoff has total occupation at least \(K+1\),

\[
\|(1-P_K)\psi_j\|\le \frac{d_j}{\sqrt{K+1}}.
\]

The compressed ordinary branch is the product \(P_KK_{s_m}\cdots P_KK_{s_1}\Omega\). Telescope this against the uncut product, using that all \(P_KK_s\) are contractions. Its vector error is at most \(B_\pi/\sqrt{K+1}\). Each ordinary-outcome probability error is at most twice that. Summing over \(2^m\) ordinary branches gives the constant in (R5).

The total probability of all histories containing overflow equals one minus the probability of the ordinary histories. Hence it is no larger than the sum of the ordinary probability errors. Extending the uncut distribution by zero on overflow histories,

\[
\operatorname{TV}(p_{\pi,\eta}^{a,K},p_{\pi,\eta}^{a,\infty})
\le \min\!\left(1,\frac{C_{\rm occ}}{\sqrt{K+1}}\right).
\tag{R6}
\]

This includes arbitrary later outcomes following an overflow. It uses no bound on the post-overflow state. The displayed constant is deliberately loose.

Thus the first P4 constant depends on the lattice spacing and size, \(L,c,\eta\), the detector smearing and strength through the displacement norms, and the history length. Phases and time gaps affect the probabilities, but this upper bound is independent of their values. The cutoff basis is fixed as above. The proof is for the vacuum preparation, not a supremum over all states.

### Spatial error

The second P4 constant is independent of \(a,K\), for odd \(N\ge3\), \(a=L/N\), and the fixed data above. Its dependence can be exposed through Fourier moments. Use an orthonormal torus Fourier convention and define

\[
\begin{aligned}
h_j(\mathbf k)&=\|\widetilde f_{E,j}(\mathbf k)\|
                 +\|\widetilde f_{B,j}(\mathbf k)\|,\\
M_{4,j}^2&=\sum_{\mathbf k\ne0}|\mathbf k|^4 h_j(\mathbf k)^2,
\quad k_{\min}=2\pi/L,\quad a_0=L/3.
\end{aligned}
\tag{R7}
\]

The smoothness assumption guarantees these moments. The frequency and covariance bounds in the bridge give a finite spectral constant \(b(c,\eta)\), independent of \(a\), such that the error in each covariance or commutator pairing is bounded by

\[
b(c,\eta)a^2(k_{\min}^{-1}+cT+a_0)M_{4,i}M_{4,j}.
\tag{R8}
\]

To see the dependencies, the amplitude error is bounded by a constant times \(a^2|\mathbf k|^3\), and the time phase error by a constant times \(a^2cT|\mathbf k|^4\). Cauchy--Schwarz and \(|\mathbf k|\ge k_{\min}\) give the first two terms of (R8). Outside the Brillouin zone, \(|\mathbf k|>\pi/a\), so the omitted covariance tail is bounded by
\((a/\pi)^3M_{4,i}M_{4,j}\le a^2(a_0/\pi^3)M_{4,i}M_{4,j}\).
The frequency interval \([(2/\pi)c|\mathbf k|,c\sqrt{1+4\eta}|\mathbf k|]\) and the stiffness estimate bound \(b\) using only \(c,\eta\). No derivative of a potentially degenerate polarization eigenvector is required. Matrix square roots on the transverse plane suffice.

Finite Weyl products have expectations \(\exp(-v/2+i\theta)\), with \(v\ge0\). Their differences are bounded by the differences of \(v/2\) and \(\theta\). Expanding the finite histories therefore gives, for example, a conservative spatial constant of the form

\[
C_{\rm sp}
 =2^{m+4}(m+1)^2\,b(c,\eta)
   (k_{\min}^{-1}+cT+a_0)
   \left(\sum_j M_{4,j}\right)^2 .
\tag{R9}
\]

The harmless enlargement of \(b\) covers the real and imaginary pairings and Fourier conventions. No optimal numerical value is claimed. Equations (R7)--(R9) state the input dependence rather than hide it in a protocol label. For a finite collection of protocols take the largest constant.

The amended P4 statement is

\[
\operatorname{TV}(p_{\pi,\eta}^{a,K},p_\pi^{\rm Maxwell})
\le \min\!\left(1,\frac{C_{\rm occ}(\pi,a,\eta,L,c)}{\sqrt{K+1}}
                       +C_{\rm sp}(\pi,\eta,L,c,T)a^2\right),
\quad
\lim_{a\to0}\lim_{K\to\infty}\operatorname{TV}=0.
\tag{R10}
\]

This is the asserted **ordered-limit theorem**. A joint sequence is licensed if
\(C_{\rm occ}(\pi,a,\eta,L,c)/\sqrt{K(a)+1}\to0\). The existence of some sufficiently fast diagonal choice follows, but the theorem makes no claim for an arbitrary simultaneous removal of the cutoffs. This review does not promote an unproved uniform occupation estimate to a theorem.

There is no asserted uniformity as \(L\) or \(T\) diverge, \(\eta\) grows, detectors approach point support, strengths or history lengths grow, preparations change without number-tail control, or conditioning probabilities approach zero. In particular this is a fixed-torus result, not the thermodynamic limit. The vacuum rule is the same for each \(\eta\), although the full vacua need not be identical.

## 3. What survives of causal commutation

Detector-probability convergence alone does not prove an exact local algebra at finite dimension. One can, however, establish a stronger, restricted statement here using the uncut Weyl relations.

For two fixed smooth probes at the specified times write
\([F_a(f),F_a(g)]=i\sigma_a(f,g)1\). The pairing estimate gives
\[
|\sigma_a-\sigma_{\rm M}|\le C_{\sigma,f,g}a^2 .
\tag{R11}
\]
Take smooth spatial smearings whose time-stamped supports are spacelike separated in the flat torus metric. This includes checking the periodic images, not just a chosen coordinate separation. Maxwell field-strength microcausality gives \(\sigma_{\rm M}=0\).
This is adopted continuum field theory. The omitted zero-mode electric and magnetic fluxes commute with all field strengths, so fixing that sector does not remove the causal commutator. A potential or Wilson-loop algebra would require separate treatment. [Dappiaggi and Lang, Section 3](https://arxiv.org/abs/1104.1374)

Before occupation truncation, the Weyl relation implies
\[
\|[W_a(f),W_a(g)]\|
  =2|\sin(\sigma_a/2)|
  \le C_{\sigma,f,g}a^2 .
\tag{R12}
\]
This norm is within each uncut lattice representation. It does not require identifying two full Hilbert spaces across different spacings.

After compression let \(P=P_{K,\eta}\), \(Q=1-P\), \(A=PW_a(f)P\) and \(B=PW_a(g)P\). The exact identity is
\[
[A,B]-P[W_a(f),W_a(g)]P
 =-PW_a(f)QW_a(g)P+PW_a(g)QW_a(f)P .
\tag{R13}
\]
Because \(P\Omega=\Omega\), it follows that
\[
\boxed{\;
\|[A,B]\Omega\|
\le \frac{\|z_f\|+\|z_g\|}{\sqrt{K+1}}
          +C_{\sigma,f,g}a^2 .
\;}
\tag{R14}
\]
Thus the compressed commutator vanishes on the prepared vacuum in the same ordered limit. At fixed \(a\), it also converges on each fixed finite coherent superposition. For a vector \(\psi\) in that class, replace the numerator by bounds for
\(\sqrt{K+1}(\|QW_f\psi\|+\|QW_g\psi\|+2\|Q\psi\|)\).
Uniformity over arbitrary cutoff-dependent boundary states is not asserted. Normalizing a rare postselected branch introduces its probability denominator.

The same reasoning applies to bounded polynomials of the ordinary compressed Kraus operators. No operator-norm convergence for the complete compressed local algebra, or its overflow square roots, follows from (R14).

### Operational no-signalling in the declared state class

Compare the finite experiment in which the first detector is applied and its outcome ignored, including overflow, with the experiment in which it is absent. Keep the second detector spacelike separated. The corresponding continuum unconditioned probabilities agree by microcausality and instrument normalization. Marginalization cannot increase total variation, so
\[
\operatorname{TV}(p_B^{a,K}\text{ with }A,\,
                  p_B^{a,K}\text{ without }A)
\le \epsilon_{AB}(a,K)+\epsilon_B(a,K),
\tag{R15}
\]
where each \(\epsilon\) is its bound from (R10). This is a causal-response limit for the declared vacuum experiments. Conditioning on a recorded outcome is a different question.

### Why an all-state claim would be wrong

Let two normal modes have canonical coordinates \(q_\pm,p_\pm\), and define two commuting spatial variables
\(q_1=(q_++q_-)/\sqrt2\), \(p_2=(p_+-p_-)/\sqrt2\).
After a product occupation cutoff at \(K\), with \(J_K=|K\rangle\langle K|\),
\[
[q_{1,K},p_{2,K}]
=\frac{i(K+1)}2(1\otimes J_K-J_K\otimes1).
\tag{R16}
\]
Its norm is \((K+1)/2\) for \(K\ge1\), despite \([q_1,p_2]=0\) before compression. This example uses unbounded fields and directly refutes automatic finite-cutoff locality. The replay additionally tests bounded commuting Weyl probes of these same two variables. Their compressed commutator is nonzero at finite \(K\), while its action on the vacuum decreases with the occupation cutoff.

The bridge now establishes joint detector probabilities, state-controlled commutator convergence, and (R15). It establishes neither exact finite-dimensional spatial commutation nor a uniform all-state causal algebra. The continuum causal structure, torus metric, and detector support map are inputs to these statements.

## 4. Decide the selection question

The declared GRUT postulates admit both \(\eta=0\) and \(\eta=1\). The uncut spatial Hamiltonians both have positivity, unitarity, gauge invariance, cubic symmetry and finite interaction range. Both admit the same axial calibration and ground-state rule, and both finite cutoff families have the fixed-probe Maxwell limit. The auxiliary occupation compression has only the guarantees established above. None of these conditions, or their conjunction, entails \(\eta=0\). This is an explicit countermodel argument within the supplied architecture.

Three possible extra rationales do not close that inference.

| Proposed rationale | What it actually requires |
|---|---|
| The simplest action | A new minimization or prior over actions, a notion of complexity, and a reason nature obeys it. A numerical choice of discretization does not supply those. |
| All observers must agree | Transporting the same physical experiment between representations already preserves its probability in both models. Exact microscopic boost symmetry at a nonzero cubic spacing would challenge both. Their common continuum limit recovers the relevant Maxwell symmetry. |
| No polarization splitting | This selects \(\eta=0\) only within the displayed one-parameter family if imposed at finite physical spacing. It is an additional empirical or theoretical requirement and does not uniquely select the plaquette action among local gauge actions in general. |

The last point can be demonstrated in a broader candidate action class using the same fields, geometry, detector dictionary and vacuum rule. This illustration is not asserted to be a member of the frozen two-term generator table. Add instead the positive cubic-invariant magnetic energy
\[
\Delta H_\zeta
=\frac{\zeta c^2a^7}{2}
  \sum_x\sum_{i<j}\sum_\ell (D_i^+D_j^+B_\ell)^2,\qquad\zeta\ge0.
\tag{R17}
\]
It multiplies both transverse squared frequencies by
\(1+\zeta a^4\sum_{i<j}\hat k_i^2\hat k_j^2\).
It is gauge invariant, vanishes on all axial modes, and adds no polarization splitting. Nevertheless it changes off-axis dispersion at finite spacing and vanishes in the same Maxwell limit. Thus even adding no birefringence would not make the existing calibration determine a unique local gauge action without an additional restriction on allowed terms. This counterterm is an analytic witness in that larger class, not a newly adopted GRUT sector or a change to the frozen comparison.

**Decision.** Retain the two Hamiltonians as an exact demonstration of nonselection under the present postulates. Do not adopt plaquette-only dynamics as a physical GRUT law. Simplicity can justify using it for computation. A claim about nature would require an independently stated selection principle or discriminating empirical input, together with identification of a physical scale. The photon construction supplies a limited microscopic bridge and a calibrated model comparison. It supplies no selection of the actual preparation, physical access, or microscopic dynamics.

## 5. Reproduction and claim scope

From the repository root, with Python, NumPy and SciPy installed, run

~~~bash
OPENBLAS_NUM_THREADS=1 python research/reality_framework/grut_integrated_candidate_20261106/verify_photon_bridge.py --output /tmp/grut_photon_bridge_review.json
~~~

The program checks complete finite axial histories including overflow, deliberately mismatched degenerate cutoffs, a noncommuting sequential-history comparison against exact uncut coherent-state algebra, compression defects for bounded commuting Weyl probes, and the spectra of (P3) and (R17). It computes displacements by their exact number-basis matrix elements and cross-checks those against a larger oscillator exponential.

The raw results are in [PHOTON_BRIDGE_RESULTS.json](PHOTON_BRIDGE_RESULTS.json). These are numerical checks of the finite identities and failure controls. The continuum and ordered-limit statements rest on the derivations above. No laboratory data or independent implementation has been substituted for those arguments.

| Finite verification | Measured result |
|---|---|
| Complete three-read axial histories, \(K=1,2,4\), vacuum and nonstationary spectators | Largest probability difference \(1.75\times10^{-15}\) |
| Deliberately rotated degenerate cutoff at \(K=1\) | Overflow difference \(0.04598493\) |
| Three noncommuting probes versus uncut coherent-state reference, \(K=24\) | Total variation \(1.90\times10^{-11}\) |
| Commuting uncut Weyl probes, compressed at \(K=1\) and \(K=12\) | Vacuum commutator norms \(0.1140\) and \(4.10\times10^{-9}\) |
| The same bounded Weyl commutator at \(K=12\) | Operator norm \(0.5919\), so the small vacuum value is not an all-state bound |

The worst-case occupation constant in (R5) is intentionally conservative. Its clipped TV bound remains one throughout the numerical sweep, including \(K=24\). The much smaller measured errors do not turn that bound into an efficient resource estimate. The replay drops only the scalar zero-point energy from each finite Hamiltonian, which cancels from every branch probability.

| Layer | Status after this clarification |
|---|---|
| General finite operational framework | Specified |
| Exact clock, record and action realization | Constructed in the frozen separate realization |
| Free flat photon detector limit | Ordered-limit derivation with explicit cutoff convention and restricted causal-response bounds, pending independent review |
| Physical selection of microscopic dynamics | Open, with explicit nonselection witnesses |
| Interacting matter and dynamical spacetime match | Open |

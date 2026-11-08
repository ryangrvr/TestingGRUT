# Eliminating states changes the detector as well as the Hamiltonian

7 November 2026 · Conditional finite-parent construction

The M1-01 continuum comparison and D1-01 reduction missed their gates. This continuation tests a derived correction scheme within the same E=12 finite interacting parent. No state above that parent cutoff is accounted for, and no old verdict changes.

The motivation is standard effective-Hamiltonian theory, not a new GRUT postulate. Integrating out omitted states can correct a truncated Hamiltonian [Cohen et al.](https://arxiv.org/abs/2110.08273). Unitary subspace transformations provide another controlled effective-Hamiltonian construction under their stated assumptions [Bravyi, DiVincenzo and Loss](https://arxiv.org/abs/1105.0675). The single-anchor isometry below is derived directly from block equations; it is neither a claimed implementation of those papers' full methods nor a perturbatively controlled Schrieffer–Wolff expansion.

## 1. Exact block elimination

For orthogonal P,Q=I-P, write the adopted finite Hamiltonian as

$$
H=\begin{pmatrix}A&B^\dagger\\ B&C\end{pmatrix}.
$$

For an eigenvalue E outside the spectrum of C, the Q equation gives

$$
q=T(E)p,\qquad T(E)=(E-C)^{-1}B,
$$

and the retained equation is

$$
H_{\rm F}(E)p=Ep,\qquad H_{\rm F}(E)=A+B^\dagger(E-C)^{-1}B.
\tag{1}
$$

The physical norm is p^dagger G(E)p, with G(E)=I+T(E)^dagger T(E). Define the isometry

$$
J(E)=\begin{pmatrix}I\\T(E)\end{pmatrix}G(E)^{-1/2},\qquad J^\dagger J=I.
\tag{2}
$$

At a chosen anchor e, the variational operators on an ordinary orthonormal retained space are

$$
H_e=J(e)^\dagger HJ(e),\qquad W_e=J(e)^\dagger WJ(e).
\tag{3}
$$

Use the lowest eigenvalue of A to initialize e, and iterate e <- lambda_min(H_e), with a fixed iteration cap and tolerance. At an exact fixed point, the lowest eigenvector v of H_e maps to a parent eigenvector: the Q residual vanishes by the defining resolvent equation, and J^dagger times the remaining P residual vanishes; the invertible P block G^(-1/2) then forces that residual to vanish too. It is necessary to verify which eigenstate was obtained; convergence to a fixed point alone does not identify the parent's ground state.

This construction need not decouple the entire retained space. Its dynamical defect is

$$
R=HJ-JH_e=(I-JJ^\dagger)HJ.
\tag{4}
$$

Even if the anchor state's residual is zero, excited superpositions can have nonzero R. There is no energy-independent exact elimination of arbitrary dynamics implied by one successful eigenvalue.

## 2. Instrument matching

The parent's ordinary operators are K_s=(I+sW)/2 and W is a contraction. Compressing with the same J gives

$$
K_{s,e}=J^\dagger K_sJ=(I+sW_e)/2.
$$

Its deficit remains positive, so the reduced process can keep an explicit absorbing overflow outcome. The original parent overflow and additional reduction loss are distinct. Replacing W_e by PWP after changing H and the state is a different experiment. Spectrum agreement cannot justify that substitution.

For the consistent construction, the missing ordinary branch vector is

$$
K_sJ\chi-JK_{s,e}\chi=\frac{s}{2}(I-JJ^\dagger)WJ\chi.
\tag{5}
$$

For any deliberately mismatched reduced W_r, its norm is instead ||(WJ-JW_r)chi||/2. Both terms, inside and outside the retained subspace, must be counted. Replacing a compressed W by a newly unitary exponential would conceal this error, not derive a correction.

## 3. A finite-parent history bound with dynamical leakage

Let chi be a reduced branch at the start of an interval of length t, and let d bound its embedded vector difference from the corresponding parent branch. Duhamel's identity gives

$$
\|e^{-iHt}J\chi-Je^{-iH_et}\chi\|
\le\int_0^t\|R e^{-iH_eu}\chi\|du.
\tag{6}
$$

There is a convenient quadrature-free upper bound. Diagonalize H_e=V diag(e_a)V^dagger, let z=V^dagger chi, and M=V^dagger R^dagger R V. Define

$$
\mathcal B_{ab}(t)=M_{ab}\,t\,e^{i(e_a-e_b)t/2}
\operatorname{sinc}\!\left[(e_a-e_b)t/2\right],
\quad \operatorname{sinc}x=\sin x/x.
$$

This is the exact positive-semidefinite integral of the squared residual on this finite space. Cauchy–Schwarz yields

$$
\int_0^t\|R e^{-iH_eu}\chi\|du\le\sqrt{t\,z^\dagger\mathcal B(t)z}.
\tag{7}
$$

Increase d by this amount, propagate chi, then increase d for each ordinary instrument by the norm in (5), or its explicit mismatched version. All parent and reduced ordinary maps are contractions. Starting with the phase-aligned preparation error d0, iteration bounds every final ordinary branch. Its probability error obeys

$$
|p_s-\widehat p_s|\le 2\|\chi_s\|d_s+d_s^2.
$$

The single absorbing-overflow outcome has probability one minus the ordinary sum in each model. Therefore

$$
\operatorname{TV}(p,\widehat p)\le
\min\!\left(1,\sum_s[2\|\chi_s\|d_s+d_s^2]\right).
\tag{8}
$$

This bound includes preparation, dynamical leakage and detector mismatch. It can be loose. Finite-precision evaluation is not an interval-arithmetic proof; positivity, residuals and agreement with the measured finite-parent discrepancy are checked explicitly. There is no bound here on the parent's omitted continuum states, different preparations, arbitrary new operations, or rare postselected probabilities.

## 4. What happens to memory

The exact retained evolution can also be written without an energy anchor. Solving the Q Schrödinger equation and substituting it gives

$$
i\dot p(t)=A p(t)+B^\dagger e^{-iCt}q(0)
-i\int_0^t B^\dagger e^{-iC(t-u)}B p(u)\,du.
\tag{9}
$$

The discarded sector returns as a memory kernel and an initial-correlation term. At a measurement, both p and q change according to the corresponding blocks of K_s. Resetting q to zero is an additional approximation. The proposed isometry is a bounded way to test whether a simpler autonomous retained description suffices in this experiment; exact block elimination by itself does not guarantee that it does.

The connection to the integrated GRUT framework is thus explicit: physical predictions belong to preparations, dynamics and instruments together. The current framework does not impose a new law that makes the discarded information irrelevant. The test below asks whether a particular derived approximation retains enough of it.

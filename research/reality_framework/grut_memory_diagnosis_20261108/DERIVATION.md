# Detector leakage, dynamical exchange, and retained memory

8 November 2026 · E2-01 · Conditional on the unchanged finite E=12 parent

E1-01 remains frozen. Its measured improvement, strong accuracy miss and both
bound misses are unchanged. The construction below diagnoses that approximation
and tests a different, memory-retaining evolution. It supplies no continuum
error estimate or new GRUT law.

## 1. The retained image and its complement

Reuse E1's isometry J, with 72 columns, without changing its anchor. Complete it
by an orthonormal complement L. In these coordinates,

\[
\psi=Jp+Lq,\quad
\widetilde H=\begin{pmatrix}A&B^\dagger\\B&C\end{pmatrix},\quad
\widetilde W=\begin{pmatrix}w&u\\v&z\end{pmatrix}.
\]

Here A=J†HJ is E1's Hamiltonian, w=J†WJ is its detector, and B=L†HJ
is its dynamical residual in complementary coordinates. P=JJ† refers to the
dressed image, not the old free-Fock projector. This is a direct-sum decomposition
of a state space, not a tensor-product system/bath decomposition.

Within any interval starting at t0, block elimination gives

\[
i\dot p(t)=Ap(t)+B^\dagger e^{-iC(t-t_0)}q(t_0)
-i\int_{t_0}^{t}B^\dagger e^{-iC(t-u)}Bp(u)\,du.
\tag{1}
\]

At an ordinary detector outcome s=±1,

\[
p_s=\tfrac12[(I+sw)p+suq],\qquad
q_s=\tfrac12[svp+(I+sz)q].
\tag{2}
\]

Even if q starts at zero, a detector generally creates q. Resetting q after a
readout removes the inhomogeneous term of the next interval and subsequent
detector return u q. Equation (1) is exact in this finite parent; evaluating it
need not be cheaper than evolving the parent.

## 2. Complete instruments, including overflow

For any contraction X, the ordinary Kraus operators are K±=(I±X)/2. Their sum
of effects is (I+X†X)/2. Introduce

\[
D_X=\frac{I-X^\dagger X}{2}\succeq0,\qquad
K_{\rm o}=V_{\rm o}\sqrt{D_X},
\tag{3}
\]

where V_o is an isometry into an absorbing copy of the input space carrying an
overflow flag. Then K+†K+ + K-†K- + Ko†Ko=I. Keeping the hidden state of the
overflow copy gives an isometric instrument dilation. Discarding it and merging
all overflow times gives exactly the single absorbing outcome used in E1.

For compression, the deficit has the explicit identity

\[
D_w=J^\dagger D_WJ+\tfrac12J^\dagger W^\dagger(I-JJ^\dagger)WJ.
\tag{4}
\]

Thus parent overflow and additional reduction loss are distinguishable positive
terms. Replacing w by exp(i times a compressed generator) is not this instrument.
No ordinary outcome is renormalized, and no negative loss is silently clipped.

## 3. Four normalized counterfactuals

Let H0=diag(A,C), V=[[0,B†],[B,0]], and W0=diag(w,z). Define

\[
H_d=H_0+dV,\quad W_i=W_0+i(\widetilde W-W_0),\qquad d,i\in\{0,1\}.
\]

The subscript i is a binary switch, not the imaginary unit. Both Hamiltonians
are Hermitian. W0 is the pinching (W+SWS)/2, S=diag(I,-I), so it is a contraction;
both detector choices have the complete instrument (3). All four cases start
from the identical embedded E1 preparation. With both switches off the process
stays in P and reproduces E1's autonomous joint scheme. With both on it is the
parent evolution, apart from the separately measured preparation discrepancy.

For the complete history distribution P_di, define signed vectors

\[
\Delta_D=P_{10}-P_{00},\quad \Delta_I=P_{01}-P_{00},\quad
\Delta_{DI}=P_{11}-P_{10}-P_{01}+P_{00}.
\tag{5}
\]

Their sum equals P11-P00 outcome by outcome, including overflow. Report every
signed vector and its half-l1 norm. These norms are not additive percentages of
the total error: interference and cancellation matter. Switching d controls
dynamical exchange; i controls detector transfer between the images. Detector-
created q can still propagate under C when d=0. Thus d=0 does not mean that all
forms of memory vanish. This factorial attribution is conditional on these
declared interventions, not a unique decomposition of physical causation.

## 4. A branch-relevant autonomous bound

Stack all ordinary, unnormalized reduced branch vectors as columns X. The parent
and reduced ordinary instruments are contractions on this stacked space. For an
interval, the exact branch-relevant propagation defect is

\[
\epsilon_D=\|(e^{-iHt}J-Je^{-iAt})X\|_F.
\]

The joint detector's stacked defect is

\[
\epsilon_I=\frac1{\sqrt2}\|(WJ-Jw)X\|_F.
\]

Repeated triangle inequalities and contraction give
delta <= delta0 + sum epsilon_D + sum epsilon_I. If Y is the final reduced stack,
Cauchy–Schwarz gives sum_s |p_s-p_hat_s| <= delta (||X_parent||F+||Y||F).
Since the parent ordinary norm is at most one, and overflow is the complement,

\[
\operatorname{TV}\le\min\{1,\delta(1+\|Y\|_F)\}.
\tag{6}
\]

All branches are retained, including those of small probability. This bound can
be smaller than E1's branchwise accumulated bound, but it uses exact parent
propagation to evaluate epsilon_D. It is a diagnostic bound, not an efficient
certificate. Its three positive budgets are bounds, not the signed attribution
in (5). No pass is inferred before evaluating it.

## 5. Normalized memory continuation and its error

Evolve the coupled (p,q) state, retaining q across every detector, with

\[
S(h)=e^{-iH_0h/2}e^{-iVh}e^{-iH_0h/2}.
\tag{7}
\]

Each factor is unitary. In an eigenbasis of A and C the first/last factors are
diagonal. If B=U Sigma R†, the middle factor rotates the paired coordinates
(R†p,U†q) with cos(Sigma h) and -i sin(Sigma h); orthogonal coordinates remain
unchanged. No singular value is discarded. The detector remains the full
contraction W with (3), so this continuation is a normalized quantum process.
Its visible p still has 72 amplitudes but its q has 1208: there is no state-space
compression claim.

For Hermitian finite matrices the two-term Strang bound is

\[
\|S(h)-e^{-i(H_0+V)h}\|\le h^3 c_S,\qquad
c_S=\frac{\|[H_0,[H_0,V]]\|}{24}
+\frac{\|[V,[H_0,V]]\|}{12}.
\tag{8}
\]

This is the real-time bound in Childs et al., *Theory of Trotter Error with
Commutator Scaling*, [arXiv:1912.08854v3](https://arxiv.org/abs/1912.08854),
equations (146)–(150), [author-hosted full text](https://yuansu.me/publication/trotter_error/trotter_error.pdf).
It follows by writing the variation-of-parameters defect, using cancellation of
its constant and linear terms, and integrating its second derivative over the
three-time simplex. Unitary conjugations preserve the commutator norms; the
two coefficients are (1/4)/6=1/24 and (1/2)/6=1/12. This is a finite-time bound,
not just a leading BCH term.

Use Frobenius norms in code as explicit upper bounds on the operator norms.
For n steps in each interval of duration t_j, telescoping unitary products gives
error <= c_S t_j^3/n^2. Lift all complete instruments to their common isometric
dilations from (3); telescoping the history circuit and then discarding registers
cannot increase trace distance. Therefore the final history TV obeys

\[
\operatorname{TV}\le\min\{1,\delta_0+c_S\sum_j t_j^3/n^2\}.
\tag{9}
\]

No multiplication by the number of outcome branches is needed. Common exact
instruments and explicit memory avoid the autonomous model's instrument/image
error; the remaining approximation is the registered propagation splitting and
preparation discrepancy. This is a different approximation, so its bound cannot
retroactively certify E1. Operator norms, matrix factorizations and probabilities
are evaluated in ordinary floating point, not interval arithmetic. The theorem
is conditional on the exact supplied finite matrices; numerical residual checks
are reported separately.

An absorbing-reset control projects q away after each nonfinal readout, puts
the removed norm into overflow, and otherwise uses exact parent propagation.
It is a different normalized instrument, not a harmless implementation choice.
It tests the consequence of erasing cross-readout memory without inventing an
unaccounted energy or probability reservoir.

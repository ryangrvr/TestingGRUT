# Exact ledger and diagnostic limits

All fields, positions and times use the pinned pilot's m_X units. Define r²=q1²+q2², j=q1 p2−q2 p1. For each positive axial displacement delta of length d or2d, the symmetric lattice coupling is respectively c_delta=4/(3d²) or−1/(12d²). Negative displacements carry the same coefficient. The force is

\[
\dot p_i=\sum_{j\sim i}c_{ij}(q_j-q_i)-\nabla U(q_i).
\]

## Local energy and flux

Share the energy of every bond equally between its endpoints:

\[
e_i={1\over2}|p_i|^2+U_i+
{1\over4}\sum_{j\sim i}c_{ij}|q_j-q_i|^2.
\]

The global sum equals the old forward-edge energy and the Fourier quadratic form. Its local allocation differs. The negative second-neighbor coefficient can make a particular local gradient allocation negative; only the global spatial quadratic form is positive. No weighted subsystem energy is silently turned into an invariant object mass.

The exact incoming energy and charge bond rates are

\[
G^E_{ij}={c_{ij}\over2}(q_j-q_i)\cdot(p_i+p_j),\qquad
G^Q_{ij}=c_{ij}(q_{1i}q_{2j}-q_{2i}q_{1j}).
\]

Both are antisymmetric under i↔j. Direct differentiation gives

\[
\dot e_i=\sum_jG^E_{ij},\qquad \dot j_i=\sum_jG^Q_{ij}.
\]

For a fixed region, outgoing flux is minus the sum of these derivatives over its cells. Interior bonds cancel. A mask changing with time would have an additional boundary-motion term. The present snapshots do not provide a time-resolved moving-core flux or an integrated escaping population.

Book the complete interaction g² phi² r²/2 in the charged Hamiltonian:

\[
e_X={|p_X|^2\over2}+e_{\nabla X}+{g^2\phi^2r^2\over2}
+{\lambda_Xr^4\over4},\qquad
e_\phi={p_\phi^2\over2}+e_{\nabla\phi}
+{\lambda_\phi\over4}(\phi^2-v^2)^2.
\]

Then e=e_X+e_phi exactly. The energy-exchange source is

\[
W_{\phi\to X}=g^2\phi\dot\phi r^2,
\quad\dot e_X=\sum_jG^X_{ij}+W_{\phi\to X},
\quad\dot e_\phi=\sum_jG^\phi_{ij}-W_{\phi\to X}.
\]

This is reversible Hamiltonian work, not a collision or decay rate. Assigning the interaction differently would redistribute the exchange term. The total energy and total outgoing flux are unchanged. The scalar distortion that binds a core is part of this ledger; so are scalar waves. Neither is subtracted as a chosen ambient constant.

## Best common rotation

For any fixed nonnegative weight w, set

\[
I=\int w r^2\,d^3x,\quad Q=\int wj\,d^3x,\quad
\bar\omega={Q\over I}.
\]

Let J(q1,q2)=(−q2,q1). Completing the square proves

\[
K_X={Q^2\over2I}
+{1\over2}\int w|p_X-\bar\omega Jq_X|^2\,d^3x.
\]

The fitted omega uniquely minimizes that squared residual when I>0. Resolving momenta radially and tangentially gives the exact nonnegative decomposition

\[
K_X={Q^2\over2I}
+{1\over2}\int w{(q_X\cdot p_X)^2\over r^2}\,d^3x
+{1\over2}\int w{(j-\bar\omega r^2)^2\over r^2}\,d^3x.
\]

The last two terms are radial motion and angular-velocity dispersion at one time, where r>0. At an exact node r=0 their directions are undefined; the implementation retains the full kinetic energy there as a separate nonnegative term. The square-completion identity is valid also at nodes. A synthetic nodal control tests this case. They are not a particle/antiparticle split. They are measured, not deleted. A fitted frequency below m_X is a necessary diagnostic for familiar isolated Q-ball branches, not a binding proof in an interacting environment.

A fixed, stationary rotating ansatz also requires p_X=omega Jq_X, p_phi=0, F_X+omega²q_X=0 and F_phi=0, where F is the full lattice force. We report the charged force residual as

\[
\mathcal R_X={\left[\int_w|F_X+\bar\omega^2q_X|^2\right]^{1/2}
\over m_X^2\left[\int_w r^2\right]^{1/2}}.
\]

A translating, excited or multi-frequency bound object can fail this stationary ansatz. A homogeneous circular condensate passes it and is not a localized bound object. The verification explicitly includes that control. Neither success nor failure alone decides boundness.

The phase-coherence diagnostic is the magnitude of the positive-current-weighted average of (q1+iq2)/r. It too is classical coherence, not a particle probability.

## Finite static coefficient and invalid endpoint approximation

On the earlier stationary circular family,

\[
D=6(1-a^2),\quad B^2=12a^2(1-a^2),\quad
{B^2\over D}=2a^2\quad (a<1).
\]

Thus the ratio has a finite limit2, and A_eff=−1.4a² is finite. The corrected statement is that the static approximation is **nonuniform and invalid** at the endpoint, not that its output diverges. Replacing the dynamical response B²/(y+k²+D) by B²/D requires|y+k²|≪D. Its domain of validity collapses as D→0. The exact finite-response factorization and independently reproduced growth rates are unchanged. This wording supersedes the corresponding loose sentences in the sealed2335653 report and mathematics note.

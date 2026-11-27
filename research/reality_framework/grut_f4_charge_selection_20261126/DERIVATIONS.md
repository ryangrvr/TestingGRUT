# Electric charge freedom under explicit symmetry breaking

## Local construction and scope

Use the Wilsonian expansion of a four-dimensional spin gauge theory. The fields and compact group are in the [frozen protocol](FROZEN_PROTOCOL.md). Use positive kinetic terms, the usual four Dirac Yukawa terms with independent generic matrices, and $V(H)=m_H^2 H^\dagger H+\lambda(H^\dagger H)^2$ with $m_H^2<0$ and $\lambda>0$. For each allowed weight vector, add
\[
\mathcal L_{\rm br}=\frac{c_B}{\Lambda^2}\mathcal O_B+
c_k\Lambda^{4-d_k}\mathcal O_k+\mathrm{h.c.},
\qquad c_B c_k\ne 0.
\]
Here $\mathcal O_B$ is the usual dimension-six $QQQL$ invariant and $\mathcal O_k$ is specified below. The coupling magnitudes are free and can be small. Use a cutoff below strong Abelian running and any thresholds for fields omitted from the EFT. This is an order-by-order perturbative local theory around the broken vacuum. The proof concerns local operators of the declared light fields. There is no common-cutoff or ultraviolet-completion assertion over the infinite family.

## 1. Solve the compact gauge and anomaly constraints

Let $(q,u,d,\ell,e,n,h)$ be the integer $U(1)_y$ weights of the left-handed $(Q,u^c,d^c,L,e^c,n^c,H)$. The four Yukawa terms use $H,H^\dagger,H^\dagger,H$ respectively. Gauge neutrality and the weak mixed anomaly require
\[
u=-q-h,\quad d=-q+h,\quad e=-\ell+h,\quad n=-\ell-h,\quad 3q+\ell=0.
\]
All equations apply per family. In relative index conventions the color mixed anomaly is $2q+u+d=0$. The mixed gravity trace is $6q+3u+3d+2\ell+e+n=0$. The cubic trace is
\[
6q^3+3u^3+3d^3+2\ell^3+e^3+n^3
=-6h^2(3q+\ell)=0.
\]
The color cubic anomaly cancels by the two fundamentals and two antifundamentals. Weak doublets have no perturbative cubic weak anomaly, and there are four per family counting color. For the specified quotient on spin manifolds, the inherited global-anomaly input is $\Omega_5^{\mathrm{Spin}}(BG)=0$ after perturbative anomaly cancellation. It is a result for the gauge group, independent of the particular allowed representation. [Davighi, Gripaios and Lohitsiri, Eq. (4.45)](https://arxiv.org/abs/1910.11277)

Write $t_3\in\{0,\pm1\}$ for color triality and $t_2\in\{0,1\}$ for weak center parity. A representation of the cover descends if $y+2t_3+3t_2=0\bmod 6$. With $h=3$ the $Q$ condition gives $q=1+6k$ and every other congruence follows from the solved equations. Thus every member, and only those members, has
\[
y_k=y_0+6k D,\quad k\in\mathbb Z,\qquad
y_0=(1,-4,2,-3,6,0,3),\quad
D=3(B-L)=(1,-1,-1,-3,3,3,0).
\]
Here $B$ is one third on quarks and $L$ is one on leptons, with opposite assignments on their left-handed conjugates. The scalar has $B=L=0$.

The cover action has exactly the declared kernel. Trivial action on $H$ first gives weak transformation $e^{-3i\alpha}\mathbf 1_2$, so the $U(1)_y$ angle is $\alpha=m\pi/3$. Trivial action on $Q$ then fixes the color transformation to $e^{2\pi i m/3}\mathbf 1_3$ for every $k$. All fields already obey the quotient condition. The quotient acts faithfully and the weights are primitive since $\gcd(h,q)=1$. Neither a change of circle normalization nor an unfaithful quotient generates this family.

## 2. Remove connected ordinary phase symmetry for every member

The generic flavor matrices have scalar commutants for the left-handed quark and lepton flavor spaces. A phase of the single Higgs and a common phase on each species therefore exhaust ordinary connected internal flavor transformations of the matter fields consistent with the kinetic and Yukawa terms, up to gauge transformations. The Yukawa equations for these phases leave three independent coordinates $(x_Q,x_L,x_H)$. Subtract $x_H y_k/3$ to set the Higgs phase to zero. The remaining directions are $B$ and $L$.

The $B-L$ current is perturbatively anomaly-free for all $k$. The color, weak and gravity traces vanish with $D$ substituted for the Abelian weight. The Abelian mixed trace $\mathrm{tr}(D y_k^2)$ also vanishes since $\mathrm{tr}(y_k^3)$ vanishes as a polynomial in the continuous interpolation $k$, and its derivative is $18\,\mathrm{tr}(D y_k^2)$. This is the accidental continuous symmetry in the renormalizable theory. The weak anomaly violates $B+L$ and preserves $B-L$.

Use a standard color, weak and Lorentz contraction for $\mathcal O_B=QQQL$. For example with $Q_i^{a j}$ of family $i$, color $a$, weak index $j$ and spinor contractions in parentheses, use $\epsilon_{abc}\epsilon_{jm}\epsilon_{ln}(Q_i^{a j}Q_{i'}^{b l})(Q_{i''}^{c m}L_{i'''}^{n})$. Distinct quark families can be chosen so the different slots give a nonzero polynomial. Its hypercharge is $3q+\ell=0$, and $(\Delta B,\Delta L)=(1,1)$. [Weinberg](https://doi.org/10.1103/PhysRevLett.43.1566)

Define $A=n^c_\alpha n^{c\alpha}$ for one flavor and $C_\mu=H^T i\sigma_2 D_\mu H$. The latter is a weak and color singlet vector of weight $2h=6$, engineering dimension 3 and $D=0$. The dimension-three Weyl contraction $A$ has $y_k=36k$ and $(\Delta B,\Delta L)=(0,-2)$. Use
\[
\mathcal O_k =
\begin{cases}
A(C_\mu C^\mu)^{3|k|},& k<0,\\
A,& k=0,\\
A[(C_\mu C^\mu)^\dagger]^{3k},& k>0.
\end{cases}
\]
The Lorentz scalar or its conjugate compensates exactly $\pm12$ units per factor. The local polynomial is gauge neutral and has $d_k=3+18|k|$.

The derivative avoids the zero contraction $H^T i\sigma_2 H$. In unitary gauge $H=(0,(v+s)/\sqrt 2)^T$ one obtains $C_\mu=i g_2(v+s)^2 W^+_\mu/(2\sqrt2)$ with $D_\mu=\partial_\mu-i g_2\tau^a W_\mu^a/2-i g_Y Y B_\mu$. The bosonic factors are unconstrained polynomial powers. Fermions occur in just one bilinear. The leading contact at nonzero $k$ can be evaluated with free transverse charged weak boson polarizations in a common transverse direction and fermion spinors with nonzero antisymmetric contraction. At first order in $c_k$ the lowest-field contact has two fermions and $6|k|$ charged vector slots. It gives a nonzero on-shell amplitude for perturbative free vector legs in a crossing channel with enough energy. For unstable weak bosons the check concerns the charged-current pole residues in the corresponding perturbative amplitudes. This supplies a nonredundancy check for the interaction. For each fixed $k$ one can choose a weak-coupling hierarchy with the relevant particle thresholds below the cutoff. With a purely Higgs vev $C_\mu$ has zero background value, so the $k\ne0$ operator does not give a Majorana mass to an electrically charged particle.

The two breaking vectors $(1,1)$ and $(0,-2)$ have determinant $-2$. In the phase coordinates after subtracting the gauge direction the first term requires $3x_Q+x_L=0$ and the second $-2x_L=0$. Thus the two independent connected ordinary phases are absent, even at the classical action level with the nonzero breaking coefficients. This result permits exact discrete and generalized symmetries outside the stated criterion. The quantum action still has an exact local gauge symmetry. An anomaly-free EFT realization of this criterion does not determine whether the ultraviolet spectrum of a theory of gravity realizes it.

## 3. Bound the first possible local $B-L$ violation

Let $d_*(k)$ be the minimum engineering dimension of a gauge-invariant, bosonic Lorentz-scalar local polynomial with net $B-L\ne0$ in the unbroken-phase fields. Neutral coefficients and curvatures are allowed, and the charge is that of an operator monomial before adding its Hermitian conjugate. Color invariance requires integer $B$ by triality. Each Weyl slot counts either $\pm1$ toward $3B$ or $\pm1$ toward $L$. Bosonic Lorentz invariance requires an even number of Weyl slots, so $3B+L$ is even. Since the difference $(3B+L)-(B-L)=2B+2L$ is even, nonzero $\Delta(B-L)$ has magnitude at least 2, or $|D_{\cal O}|\ge6$. [Related operator parity in canonical SMEFT](https://doi.org/10.1016/j.physletb.2016.05.050)

For $s=1/2$ and $3/2$, each elementary field and its conjugate obey
\[
|y_0-sD|\le 3 d_{\rm eng}.
\]
The finite check is on the [input table](PARAMETERS.json). Derivatives and gauge curvature are neutral in this inequality and only increase dimension. Add the inequalities over the monomial and use the exact equation $y_0(\mathcal O)=-6kD_{\cal O}$. Then
\[
d_{\rm eng}(\mathcal O)\ \ge\ \frac{|6k+s|\,|D_{\cal O}|}{3}\ \ge\ 2|6k+s|.
\]
Use $s=3/2$ at $k\ge0$ and $s=1/2$ for $k<0$. With the explicit upper witness this yields
\[
\begin{array}{ll}
12k+3\le d_*(k)\le18k+3,&k\ge0,\\
12|k|-1\le d_*(k)\le18|k|+3,&k<0 .
\end{array}
\]
At $k=0$ both sides equal 3. At $k\ne0$ these are bounds. At every finite dimension cap there are full EFTs without exact connected ordinary phases whose truncations still preserve continuous $B-L$. This separates an exact symmetry premise from the stronger premise that some low-dimensional operator must violate it. Within this class, an independently established operator of dimension at most 10 with nonzero $B-L$ would imply $k=0$ since the smallest noncanonical lower bound is 11. In particular, a gauge-invariant dimension-five Weinberg term has weight $2(\ell+h)=-36k$ and therefore forces $k=0$ within this class. The Majorana bilinear forces the same result. No such term is used to select $k$ in the test.

## 4. Physical equivalence, radiative stability and running

The vev fixes the unbroken generator to $Q_{\rm em}=T_3+y_k/6$, with unit positive weak boson charge. The electric flux and soft electromagnetic couplings of identifiable states measure representation ratios. In particular
\[
(q_u,q_d,q_e,q_\nu)=\left(\tfrac23+k,-\tfrac13+k,-1-3k,-3k\right).
\]
Additivity in each strong-interaction charge sector gives $q_p=1+3k$ for $uud$ and $q_n=3k$ for $udd$. The ratios differ physically, including $q_n/|q_e|=3k/|1+3k|$. Hydrogen as proton plus electron is neutral across the family. The quantitative comparison with a measured Coulomb charge is in the [reality test](REALITY_TEST.md).

Gauge-compatible local changes of field variables and renormalization preserve the action of the compact gauge group and electromagnetic Ward identities. The generator normalization drops out of charge ratios. The operator bound concerns the first gauge-invariant selection-rule violation in the canonical low-energy field and derivative expansion. Counterterms with the same light fields obey the same charge and spin rules. Wilson coefficients and basis choices can depend on the scheme. The proof does not bound anomalous dimensions in a strongly coupled fixed point, rates in the broken theory after heavy thresholds are removed, or interactions with extra charged UV fields. Setting every $B-L$ breaking coupling to zero in the local EFT restores the accidental anomaly-free phase, so perturbative mixing preserves that limit. No numerical lower bound on the size of breaking follows.

As a conventional check use $Y=y_k/6$ and three active families, plus one active Higgs doublet. The ordinary one-loop Abelian coefficient for canonical weak coupling is
\[
\mu\frac{dg_Y}{d\mu}=\frac{g_Y^3}{16\pi^2}b_Y+\cdots,\qquad
b_Y=\frac23\sum_{\rm Weyl}Y_i^2+\frac13\sum_{\rm complex\ scalar}Y_i^2
=\frac{41+192k+576k^2}{6}.
\]
This is the leading dimension-four running in the weak-coupling region above the relevant masses and below the cutoff. A different generator normalization rescales $g_Y$ and $b_Y$ together. RG consistency evolves the coupling supplied at a matching scale. It leaves the representation label and matching value independent. In particular the EFT construction can lower its Abelian coupling for large $|k|$ to retain a perturbative range. A shared ultraviolet embedding, maximum charge, Landau-pole margin at fixed coupling or minimum breaking rate would add premises.

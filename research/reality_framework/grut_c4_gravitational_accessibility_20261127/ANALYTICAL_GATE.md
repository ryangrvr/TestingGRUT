# C4 analytical gate

The candidate is the unchanged original Einstein–scalar action at the pinned C3 benchmark. Use the C2 soliton at exact charge $Q=39{,}944$ as the leading classical source. The C3 branch is the base and the C2 enclosure is an input. The calculation derives leading metric scales and one bound on static scalar rearrangement. The derivative and Newton bounds in C2 are for the leading flat scalar theory. Solving the coupled Einstein–scalar background, certifying gravitating stability and calculating a fully matched relativistic tidal response would require additional work. [Action](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md) · [C2 certificate](../grut_c2_continuum_stability_20261127/REPORT.md) · [Protocol](PROTOCOL.json)

## Conserved source and spherical metric

In $\rho=g_Xvr$ and $\tau=g_Xvt$, take $X=vF e^{-i\Omega\tau}/\sqrt{2}$, $\phi=vS$ and $D=0$. Orthogonal-frame stress components divided by $g_X^2v^4$ are
\[
K_t=\tfrac12\Omega^2F^2,\quad K_r=\tfrac12(F'^2+S'^2),
\quad W=\tfrac14(S^2-1)^2+\tfrac12 S^2F^2+\tfrac1{40}F^4,
\]
\[
e=K_t+K_r+W,\quad p_r=K_t+K_r-W,\quad p_t=K_t-K_r-W.
\]
The time-radial flux vanishes. The flat field equations give $p_r'+2(p_r-p_t)/\rho=0$, so $\int_0^\infty\rho^2(p_r+2p_t)d\rho=[\rho^3p_r]_0^\infty=0$. As a result, the active source $J=e+p_r+2p_t=4K_t-2W$ has the same integrated monopole as $e$. Pressure affects the nonuniform interior lapse. With $I(\rho)=\int_0^\rho s^2e(s)ds$ and $\alpha=v^2/(2M_P^2)$, the leading areal metric has $2Gm/(rc^2)=2\alpha I(\rho)/\rho$ and
\[
\nu'(\rho)=\alpha\left[I(\rho)/\rho^2+\rho p_r(\rho)\right],\qquad \nu(\infty)=0.
\]
Here $g_{tt}=-e^{2\nu}$. The exact spherical constraint would divide the right-hand side by $1-2Gm/(rc^2)$ and use the curved-saddle stress. $M_\infty c^2=(v/g_X)4\pi I(\infty)$ at leading classical order. The C1 energy comparison includes the scalar one-loop correction to flat rest energy. It supplies neither a renormalized local one-loop stress nor the gravitational binding correction to an exact ADM mass.

Birkhoff's theorem applies on spherically symmetric vacuum regions. With mass and the normalization of exterior time matched, static exterior vacuum observables are Schwarzschild observables. The exact leading scalar profile is noncompact. A vacuum profile on an open exterior interval would instead give the vacuum solution everywhere by uniqueness of the regular radial equations away from the origin. Use the decaying actual source to compare the asymptotic exterior. [Vacuum theorem](https://arxiv.org/abs/0908.4110) · [Original weak gravity](../grut_candidate_1_20261018/WORKED_REGIME.md)

## Exact-profile mass and scalar-tail envelope

Put $y=(u,w)=(\rho F,\rho(S-1))$. The C2 radius $\eta$ bounds $\|y-\hat y\|_2$ and $\|\smash{y''-\hat y''}\|_2/16$. It gives $\|(y-\hat y)'\|_2\le4\eta$ and $\|y-\hat y\|_\infty\le3\eta$. In weighted integrals the gradient energy can use $\int \rho^2(F'^2+S'^2)=\|y'\|_2^2$ by integration by parts. Exact-profile mass is enclosed using this estimate, the frequency interval and the interval integral of the predictor potential. The potential perturbation obeys
\[
\left|\int \rho^2\big(W(F,S)-W(\hat F,\hat S)\big)d\rho\right|
\le 10(\|\hat y\|_2+\eta)\eta.
\]
This follows from $|\nabla W|\le 10\sqrt{F^2+(S-1)^2}$ on the full interpolating tube. The two gradient components satisfy $|W_F|\le 4.9|F|$ and $|W_S|\le6(|F|+|S-1|)$ there.

At the C2 taper radius $R=60.125$ the exact half-line bounds give $|y|\le a e^{-\kappa(\rho-R)}$ and $|y'|\le b e^{-\kappa(\rho-R)}$, with $\kappa=1/3$. Since $|p_r|,|p_t|\le e$, an envelope for stress and exterior mass is
\[
\rho^2 e(\rho)\le B e^{-2\kappa(\rho-R)},\quad
B=\tfrac12(b+a/R)^2+3a^2,\qquad
I(\infty)-I(\rho)\le \frac{B}{2\kappa} e^{-2\kappa(\rho-R)} .
\]
The gradient part uses $y'-y/\rho$. The tube caps bound the kinetic and potential terms by $3|y|^2$ after multiplication by $\rho^2$. Comparing the areal lapse to the matched static vacuum field at the same order, the tail bounds
\[
|\delta\nu(\rho)|\le
\frac{\alpha B}{2\kappa}\left(\frac{1}{2\kappa\rho^2}+\frac1\rho\right)e^{-2\kappa(\rho-R)}.
\]
A similar bound on $|\delta\nu'|$ follows by omitting the final factor $1/(2\kappa)$ outside the parentheses. The bounds concern the leading classical source and the asymptotically normalized lapse.

The conserved stress yields a sharper check. It gives $\int_0^\rho x^2J\,dx=I(\rho)+\rho^3p_r(\rho)$. Thus the same lapse obeys the radial Poisson equation for $J$. With $|J|\le4e$,
\[
|\delta\nu(\rho)|=\left|\alpha\int_\rho^\infty x^2J(x)(1/\rho-1/x)\,dx\right|
\le\frac{4\alpha B}{(2\kappa)^2\rho^2}e^{-2\kappa(\rho-R)} .
\]
The corresponding gradient bound is $4\alpha B e^{-2\kappa(\rho-R)}/(2\kappa\rho^2)$.

## Sole dynamical screen

Screen the static scalar rearrangement under a slowly varying external quadrupolar tide. Use a real harmonic normalized by $\int Y_{2m}^2d\Omega=1$ and impose $\Phi_{\rm ext}=t\,\rho^2Y_{2m}\chi(\rho)$ in the target region, with $0\le\chi\le1$ and $\chi=1$ across the occupied region. The window restricts the polynomial tide to a weak near-zone field. In the leading weak metric $g_{00}=-(1+2\Phi_{\rm ext})$, $g_{ij}=(1-2\Phi_{\rm ext})\delta_{ij}$, the spatial-gradient prefactor is unchanged at linear order. Lapse and volume change the stationary Routhian by $(v/g_X)\int \Phi_{\rm ext}J\,d^3\rho$. The first-order total charge change integrates to zero for the quadrupole. The perturbation is therefore governed by the C2 $\ell=2$ amplitude form even after the charge constraint is imposed.

Write $b_2(\rho)=\rho^3\chi(\rho)\nabla_{F,S}J$, and let $\mathcal I_2=\int \rho^4Y_{2m}\chi J\,d\rho d\Omega$ be the dimensionless conjugate energy moment. Static matter rearrangement gives
\[
-\partial_t\mathcal I_2\big|_{\rm rearr.}=\langle b_2,L_{+,2}^{-1}b_2\rangle .
\]
The symbol $t$ here is a tidal amplitude. The C2 amplitude form $L_{+,1}$ is nonnegative. Its matrix potential has norm at most 36 by the field tube caps, and its exterior lower form is at least $1/8$ in the units of the static operator. With $\theta=2/(36R^2)$, combine $\theta L_{+,1}$, $(1-\theta)L_{+,1}\ge0$ and the additional $4/\rho^2$. Inside $R$ the lower bound is $2/R^2$. Outside it is $\theta/8$. This gives
\[
\lambda_{*,2}=\frac{1}{144R^2}\le L_{+,2},\qquad
|\partial_t\mathcal I_2|_{\rm rearr.}\le\frac{\|b_2\|^2}{\lambda_{*,2}}.
\]
The same conservative stiffness bound applies to the phase $\ell=2$ form using the exact phase ground state and centrifugal ordering. It does not specify internal resonances. The force-norm enclosure uses C2 Bernstein profile hulls expanded by $\min(12\eta,3\eta/\rho_{\rm left})$ in each cell. In the exterior, $|\nabla J|^2\le484|y|^2/\rho^2$, so the remaining integral is bounded by $484a^2\int_R^\infty\rho^4e^{-2\kappa(\rho-R)}d\rho$. This permits an interval bound without response inversion.

The bound controls the specified material contribution only. A relativistic Love number also requires metric response, contact terms and gauge-invariant multipole matching with the external source. Established theory permits such tidal differences between matched-mass objects with different internal dynamics. The benchmark's conventional quartic FLS realization with the same matched stress and scalar equations has the same response. This is the matched-comparator obstruction relevant to the feasibility decision. [Tidal response convention](https://arxiv.org/abs/0906.1366) · [Conventional scalar and compact-object comparison](https://arxiv.org/abs/1701.01116)

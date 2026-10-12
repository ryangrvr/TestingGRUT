# SCOUT-RG01 — relational matter gaps, extensive equilibrium, and an acoustic causal comparator

**2026-10-11.** A new physical hypothesis tested independently after the [fixed-bond quantum entropy obstruction](https://github.com/ryangrvr/TestingGRUT/blob/2e2c9e98005eb9d42990b9528ee2a713001217e0/research/reality_framework/cp11_fixed_bond_entropy_no_go/REPORT.md). This research packet is on a new branch forked from 2e2c9e9. Frozen CP-11, GRUT-RAI and Stage-3 Card 3 remain untouched. **This is a constructive 1D counterexample to the need for independent graph qubits, not a derived 3+1-dimensional universe or original physics principle.**

## Question

CP-11 treats all possible pairwise edges as physically independent configurations. A fixed-bond inventory can enforce sparsity but does not defeat the O(N log N) combinatorial surplus of nonspatial graphs if energy differences remain O(N). What happens if **connections are not microscopic Hilbert-space degrees of freedom at all**, but relational consequences of a smaller matter configuration space?

## Law and what we supply

Take N identical impenetrable particles with a compact **one-dimensional cyclic phase-coordinate topology S1**; the angular circumference L is variable. Order particles cyclically and denote their N successive physical gaps by \ell_i>0, with L=\sum_i\ell_i. Declare
\[
\boxed{H_{\rm fixed\ L}=\sum_i p_i^2/(2m)+g\sum_i1/\ell_i+PL,\qquad g,P,m>0.}
\]
For the fixed-L acoustic dynamics, PL is a constant. The pressure P is an externally supplied physical tension, not an Einstein equation or cosmological constant. The gap interaction is *nearest neighbor in the supplied cyclic phase topology*. It cannot be presented as emergent three-dimensional locality. A fully autonomous quantum operator for varying L, with its kinetic and domain/boundary conditions, is not implemented; volume minimization and the thermal isobaric ensemble are separate well-defined benchmarks. Full quantum results are limited to the harmonic fluctuation theory.

**No graph occupation qubits exist.** For given matter phases, neighbor connectivity is a derived successor relation. Identical-particle permutation labels are gauge for this benchmark. The dimension/topology of the phase coordinate, the exclusion of particle crossings, the repulsive potential, clock time, Hilbert-space structure and tension are OPENLY SUPPLIED.

## Result 1 — stable physical size without a separately chosen lattice spacing

The elementary convex inequality (Jensen/Cauchy-Schwarz) gives, at given L,
\[
\sum_i1/\ell_i\ge N^2/L
\]
with equality only when every \ell_i=L/N. The effective classical minimum is
\[
E(L)=gN^2/L+PL.
\]
Minimizing it over L yields
\[
\boxed{L_*=N\sqrt{g/P},\quad
a_*=\sqrt{g/P},\quad E_{\min}=2N\sqrt{gP}.}
\]
Thus the circumference grows linearly with N, the local gap is finite and N-independent, and the minimum energy is extensive. This is **not** spontaneous selection of the circle's topology: S1 is a property of the original particle coordinate.

The early counterfactual of holding L=2pi fixed while increasing N had a \Theta(N²) harmonic zero-point energy. It was a bad thermodynamic limit. The variable-L density law above repairs that particular scaling without tuning m, g, P or hbar with N.

## Result 2 — the same gap energy produces a wave

Put x_i=i a_*+u_i, at fixed L=L*, and retain quadratic terms:
\[
V=V_{\min}+\frac{g}{a_*^3}\sum_i(u_{i+1}-u_i)^2+O(u^3).
\]
The Hessian is \(K=(2g/a_*^3)L_{\rm cycle}\), derived independently by finite differencing the microscopic inverse-gap force. Its Fourier eigenfrequencies are
\[
\boxed{\omega_k^2=\frac{8g}{m a_*^3}\sin^2(\pi k/N),\qquad
q_k=\frac{2\pi k}{N a_*}.}
\]
Therefore for k/N->0,
\[
\boxed{\omega_k^2=c^2q_k^2+O(q_k^4a_*^4),\quad
c=\sqrt{\frac{2g}{m a_*}}.}
\]
The source/detector coupling range in the harmonic Hamiltonian is precisely the induced cyclic-neighbor relation. The group velocity is bounded by c in its large-N band; this is an **acoustic, approximately local, nonrelativistic** wave, not a derived universal Einstein light cone.

There is a sharper EXACT impulse test. A unit momentum impulse at particle 0 gives, for an infinite uniform ring/chain in the harmonic limit, the velocity response at relational hop distance r,
\[
\boxed{m\,\partial_t u_r(t)=J_{2r}(\Omega t),\qquad
\Omega=\sqrt{8g/(m a_*^3)}.}
\]
The finite-N ring response is its periodic-image sum. The Bessel function \(J_{2r}(z)\) is strongly suppressed when 2r is large relative to z, consistent with a front at r a_*≈ct; there are no perfectly sharp relativistic signal cones at finite lattice spacing. Independently computing the full normal-mode response at N=128 agrees with J_{2r} at all 15 tested (r,t) pairs to <1e-11.

The adversarial control is a permutation-symmetric **all-pairs** logarithmic phase repulsion on a *supplied fixed circle*. It also favors regular angular separation, but its Hessian couples any two sites immediately: the momentum-kick response has a nonzero offsite t³ term for every separation, whereas the nearest-gap Hamiltonian first couples hop distance r at order t^(2r+1). Its long-wave oscillator spectrum at that fixed circle scales as \omega_k²∝k(N-k)/N, giving \omega∝sqrt(k) rather than k for k≪N. This comparator does **not** keep the same thermodynamic pressure law; it isolates the role of microscopic interaction locality, not a universal superiority of this specific inverse-gap potential.

### Independent propagation-derived spectral dimension

Using the same induced-neighbor Laplacian, define the heat-return probability
\[
P_N(s)=\frac1N\sum_{k=0}^{N-1}\exp[-4s\sin^2(\pi k/N)],
\qquad d_s(s)=-2\,\frac{d\ln P_N}{d\ln s}.
\]
For N=256, evaluated at diffusion time s=5, 20, 80 (well before finite-ring saturation), the numerical spectral dimensions are approximately **1.0280, 1.0064, 1.0016**, respectively. Thus the transport operator correctly reconstructs the same **one-dimensional** local geometry selected as a cyclic phase ordering. This is an **operational correspondence test**, not a derived spatial dimension: the underlying S1 cyclic topology and nearest-phase-neighbor interaction were already supplied.

## Result 3 — harmonic quantum energy is extensive

At fixed equilibrium a*, quantizing the N−1 nonzero harmonic normal modes gives
\[
\boxed{
E_{\rm zpe}(N)=\frac{\hbar\Omega}{2}\cot\left(\frac{\pi}{2N}\right)
=\frac{\hbar\Omega}{\pi}N+O(N^{-1}).
}
\]
Local neighbor-gap fluctuation variance is
\[
\boxed{\langle(u_{i+1}-u_i)^2\rangle
=\frac{2\hbar}{mN\Omega}\cot\left(\frac{\pi}{2N}\right)
\longrightarrow\frac{4\hbar}{m\pi\Omega}.}
\]
Unlike the fixed-circumference trial, these local quantum fluctuations and the zero-point energy **do not diverge per particle** as N grows. For the declared control g=P=hbar=1, m=100: a*=1; c=0.1414213562; \Omega=0.2828427125; \(E_{\rm zpe}/N\to0.0900316316\), with neighbor-gap rms \(\to0.2122\). The small-fluctuation harmonic approximation is reasonably controlled at these parameters, but is not proof of the global nonlinear quantum ground state.

The zero-point term also shifts the equilibrium radius. The self-consistent **harmonic/one-loop static** energy density is
\[
e_{\rm harm}(a)=g/a+Pa+\frac{\hbar}{\pi}\sqrt{\frac{8g}{ma^3}},
\]
whose minimum for these controls is \(a_{\rm harm}\approx1.0634641\), \(e_{\rm harm}\approx2.0858812\), rather than a=1 and e=2. These are approximate effective static quantities and **not** the spectrum of a constructed, fully dynamical quantized length operator. Infrared positional fluctuations in 1D may grow with size even while local gap fluctuations stay finite.

## Result 4 — finite positive temperature has a finite gap distribution

Gauge fix the overall cyclic translation and use the declared classical isothermal-isobaric configurational measure \(\prod_i d\ell_i\). The gap distribution is then exactly
\[
p_\beta(\ell)=Z_1^{-1}\exp[-\beta(g/\ell+P\ell)],\qquad
\boxed{Z_1=2\sqrt{g/P}\,K_1(2\beta\sqrt{gP})},
\]
where K_1 is the modified Bessel function. Consequently,
\[
\boxed{\langle\ell\rangle_\beta=
\sqrt{g/P}\,\frac{K_2(2\beta\sqrt{gP})}
{K_1(2\beta\sqrt{gP})}.}
\]
Both positive and inverse gap moments are finite for every fixed beta>0 and g,P>0. Thus \(\langle L\rangle=N\langle\ell\rangle\) grows extensively in the *classical isobaric model*. Numeric integration of Z_1 and independent Bessel evaluation agree to <1e-10. With g=P=1, mean gap is 2.699484 (beta=.5), 1.814308 (beta=1), 1.258980 (beta=3), 1.075893 (beta=10).

This does not imply an ordered crystal persists at positive T: cumulative position fluctuations in 1D can be unbounded, and the quantum finite-T ensemble has not been solved. It demonstrates a well-defined local gap fluid, **rather than selecting a spatial ring from competing graph topologies**.

## Result 5 — actual configuration-space size changes

If the compact phase is discretized into q=2N available phase bins and the N particles are identical and hard-core, the number of configurations is
\[
\boxed{\binom{2N}{N}\sim 4^N/\sqrt{\pi N}=\exp(O(N)).}
\]
Each configuration induces cyclic predecessor/successor relations; independent connections are never enumerated. The unrestricted K=N simple-graph space would instead have \(\binom{N(N-1)/2}{N}\), with \(\log\sim N\log N\). This contrasts the *kinematical resources*, not two equivalent statistical ensembles: choosing S1, bin resolution q∝N, hard-core particles and permutation-gauge identification is the entire reason the graph entropy obstruction no longer applies. **A circle has been supplied at a deeper level, not selected from all other topologies.**

## Executable results

[check.py](check.py) tests N=8,16,32,64,128: variational energies, independent finite-difference force-Hessian residual (<1e-8), exact Fourier modes and finite-size errors (<1.1e-4 for the N=128 k=1 acoustic comparator), finite quantum zero-point sum and local variance, and 15 independent Bessel versus full-mode impulse responses, plus three propagation-derived spectral-dimension checks. It checks the classical isobaric one-gap integral at four temperatures against K_1 and the phase-space versus unrestricted graph state count. It also checks the one-loop static equilibrium shift. Equivalent independent standalone calculations passed before committing. The exact committed script should be run for repository confirmation; its output file results.json is generated locally.

This is a new **model-specific constructive demonstration** assembled from known 1D mechanics, 1D phonon chains and exact nearest-neighbor fluid thermodynamics. No claim of literature priority or new empirical prediction is made.

Relevant prior art:
- [Phonons of a 1D chain (standard nearest-neighbor phonon derivation)](https://materia.fisica.unimi.it/manini/dida/1D-phonons.notes.html).
- [Schulz, Wigner Crystal in One Dimension (1993)](https://doi.org/10.1103/PhysRevLett.71.1864).
- [Pustilnik & Matveev, Solitons in a one-dimensional Wigner crystal (2015)](https://doi.org/10.1103/PhysRevB.91.165416).
- [Santos, one-dimensional exact nearest-neighbor isothermal-isobaric fluids](https://ndl.ethernet.edu.et/bitstream/123456789/56503/1/Andr%C3%A9s%20Santos.pdf).

## Physical failure that remains and a specific next test

The construction **passes**: induced graph adjacency without edge qubits; a dynamically preferred average separation and circumference; extensive classical and harmonic quantum energies; a gapless bounded-velocity acoustic branch tied to the induced neighbors; and a finite positive-temperature classical local-gap state.

It **fails as a universe** because the *one-dimensional cyclic phase manifold and nearest-phase-neighbor locality were chosen*. Spatial dimension, topology, Einstein gravity, Lorentz invariance, spin-2 excitations, quantum outcome selection, clocks and identification of internal phase with physical space remain unearned. Even the fully dynamical circumference quantum Hamiltonian is still absent. Ordinary 1D phonons are not evidence for gravity.

The next physically decisive construction should not generalize the model by **assuming** S1×S1×S1 or an ambient R³ and calling it dimension selection. It must instead allow competing phase/order topologies or purely operator-algebraic relational neighborhoods under one microscopic state space, then show that a finite-dimensional locally propagating sector dominates and that its matter excitations are consistent with the same derived geometry. That is a genuinely open problem, not an arithmetic continuation of this model.

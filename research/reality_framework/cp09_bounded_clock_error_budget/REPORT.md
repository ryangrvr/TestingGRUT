# CP-09 bounded clock: matched finite-width comparator, quantum dispersion, and resource scaling

Date: 2026-10-11. Independent, additive research note. Frozen CP-09, its [moving translation-clock control](https://github.com/ryangrvr/TestingGRUT/blob/401de9b1b37a1784c25cc3485130ef5a47fafdfa/research/reality_framework/cp09_moving_clock_apparatus/REPORT.md), and its [bounded-below clock simulation](https://github.com/ryangrvr/TestingGRUT/blob/f6288d647a7da90995b6921c6c0b49ad8ce8bdbf/research/reality_framework/cp09_bounded_clock/REPORT.md) are UNCHANGED. This review forks from f6288d6, opens no CP-10/Stage-3 Card 3, and proposes no new microscopic causal law. The published numerical bounded-clock result is *imported as reported*, not independently rerun on a grid.

## 1. Separate the finite-width baseline from the bounded-clock residual

All numbers here refer to CP-09 at **N=2**, the same disjoint minimal/maximal-role source B and detector A, source-to-detector separation tau=0.4, J=1, Omega=1/4 and gamma=pi/4.

The earlier unbounded-below translation-clock Hamiltonian \(H_c=P\) implements **finite square-window** pulses exactly by its method of characteristics. Its source and pointer windows have equal width delta and integrated strengths epsilon=zeta. The later positive \(H_c=P^2/(2m)\) experiment used delta=0.15 and **epsilon=zeta=0.001**, not 0.0001 as in the earlier translation-clock report. We therefore independently rebuilt CP-09's 12-by-12 system Hamiltonian and the two 24-by-24 exact chiral pulse propagators, and **used the same epsilon=zeta=0.001** for an apples-to-apples finite-window comparison.

For initial joint system ground state and |0> pointer, the sign-difference estimator is
\[
\widehat R_\delta =
-\frac{\langle X_r\rangle_{+\epsilon}-\langle X_r\rangle_{-\epsilon}}
{4\epsilon\zeta}.
\]
The translation-clock result for two nonoverlapping square windows is computed by
\[
U_s^\pm=e^{-i(H_s\mp(\epsilon/\delta)B)\delta},\quad
U_0=e^{-iH_s(\tau-\delta)},\quad
U_d=e^{-i(H_s\otimes I_r-(\zeta/\delta)A\otimes Y_r)\delta}
\]
(with the appropriate pointer identities). This is exact for the **specified ideal translating clock**, not the dispersive clock.

| Comparison | Response |
|---|---:|
| Instantaneous ideal Kubo, R(0.4) | 0.052588513538745 |
| Translation-clock square windows delta=0.2, strengths 0.001 | 0.052159929639416 |
| Translation-clock delta=0.15, strengths 0.001 | **0.052349192601293** |
| Translation-clock delta=0.1, strengths 0.001 | 0.052482689981056 |
| Translation-clock delta=0.05, strengths 0.001 | 0.052562104906958 |
| Translation-clock delta=0.025, strengths 0.001 | 0.052581877886971 |
| Previously reported bounded clock, m=320, grid=2048, dt=.00125, delta=.15 | **0.0539630583** |

Hence the observed difference from the ideal splits *arithmetically* into
\[
\boxed{
R_{\rm bounded}-R_{\rm ideal}
=\underbrace{(-0.000239320937452)}_{\text{finite-window translation benchmark}}
+\underbrace{(+0.001613865698707)}_{\text{bounded result minus matched benchmark}}
=+0.001374544761255.
}
\]

This is a useful CONTROL COMPARISON, **not** a rigorous decomposition of independent physical error sources: the latter term combines finite clock mass, wavepacket dispersion, clock back-action, different clock transport law, periodic circle boundary, and residual lattice/timestep error. The reported single refinement at m=320 changes R by 8.98e-6 but does not certify the full discretization error bound. Finite perturbative source strengths also differ from a strict zero-strength derivative, though now matched exactly between the two apparatus calculations. The opposite signs mean the quoted 2.61% *net* response error benefits from cancellation.

A narrower pulse improves the translation-clock limit, but *cannot* be assumed to improve the positive quadratic clock at a fixed mass: a new two-dimensional (mass, window width) convergence study would be needed.

## 2. Exact freely propagating Gaussian clock resource identity

For a free clock on the line with \(H_c=P^2/(2m)\), \([q,P]=i\), set \(\hbar=1\). Its Heisenberg position is \(q(t)=q(0)+tP/m\), and
\[
\boxed{[q(t_s),q(t_d)]=i(t_d-t_s)/m.}
\]
The finite-mass clock position at separate times is thus not a pair of commuting classical clock readings.

For an initially *unchirped minimum-uncertainty Gaussian* with initial position standard deviation sigma, mean momentum p0=m (speed 1) and momentum standard deviation 1/(2sigma):
\[
\boxed{\sigma_q^2(t)=\sigma^2+\frac{t^2}{4m^2\sigma^2},\quad
\langle H_c\rangle=\frac{m}{2}+\frac{1}{8m\sigma^2},\quad
(\Delta H_c)^2=\frac{1}{4\sigma^2}+\frac{1}{32m^2\sigma^4}.}
\]
The apparatus starts near q=-0.7, so free classical travel to the source and detector centers (q=0, q=0.4) takes t_s=.7 and t_d=1.1. With sigma=.08 as in the bounded-clock code:

| mass m | sigma_q at source | sigma_q at detector | mean free kinetic energy |
|---:|---:|---:|---:|
| 20 | 0.232920 | 0.352936 | 10.976563 |
| 80 | 0.096906 | 0.117411 | 40.244141 |
| 320 | 0.081160 | 0.082835 | 160.061035 |

These are **exact infinite-line, free-packet benchmarks**, not numerical outputs from the interacting periodic-grid apparatus. Before physical interaction, the packet at m=320 is approximately narrower than the delta=.15 square window but not asymptotically sharp compared with it. At m=20 it is substantially broader. The reported full initial energies add CP-09's ground energy and small initial interaction tails to these free values.

Minimizing sigma_q²(T) over *unchirped* initial Gaussian sigma gives
\[
\boxed{\min_\sigma \sigma_q(T)=\sqrt{T/m},}
\]
so holding free arrival resolution smaller than target width delta across later flight T requires \(m\gtrsim T/\delta^2\), meaning leading free clock kinetic energy m/2 grows at least as \(O(\delta^{-2})\) in this preparation family. **This is not a universal bound for all clock states**: initially correlated/chirped packets can focus near a selected time, but maintaining accuracy at multiple source/detector readings and paying preparation/energy costs remains to be separately modeled.

## 3. A further exact finite-mass correction: clock window noncommutation

Choose two bounded Fourier-mode window probes \(W_s=\exp[i k q(t_s)]\), \(W_d=\exp[i l q(t_d)]\). The c-number Heisenberg commutator gives an exact Weyl relation
\[
W_s W_d=e^{-ikl(t_d-t_s)/m}W_d W_s,\quad
\boxed{\|[W_s,W_d]\|=2\left|\sin\frac{kl(t_d-t_s)}{2m}\right|.}
\]
The clock therefore cannot in general be represented by one nonbackreacting classical random arrival time across both windows. This is a precise **operator** diagnostic independent of the prepared packet.

For Fourier scales \(k=l=1/\delta\), define eta=tau/(m delta²) for tau=.4, delta=.15. Its values are 0.888889, 0.222222, 0.055556 for m=20, 80, 320. The corresponding Weyl commutator norms are 0.859913, 0.221765, 0.055548. These are illustrative gate-frequency scales; the **square** window has arbitrarily high Fourier harmonics. None of the norms is a bound on the measured pointer-response error and none is proof that these particular noncommuting probes were separately measured in the published apparatus.

Increasing m drives this noncommutation to zero but raises the free kinetic energy at fixed speed and momentum-space preparation. Narrowing delta increases eta unless m is increased accordingly. This exposes an independent, physically grounded finite-clock issue, beyond numerical grid convergence or the width of the square pulse.

## 4. Reproducibility and scientific boundary

[check.py](check.py) reconstructs the CP-09 N=2 Hamiltonian without importing its original modules, computes the exact finite-window translation-clock pointer expectation at the **same 0.001 source/pointer strengths** as the bounded apparatus, reproduces the ideal Kubo number, and evaluates free Gaussian dispersion/energy and Weyl window commutator scales. It was run locally before committing and passes with NumPy+SciPy.

It does NOT rerun the full positive-clock grid simulation (mass, grid, timestep sequence), verify exact continuum-grid convergence, or physically separate all the causes of the bounded-clock discrepancy. The \`R_{\rm bounded}\` value is taken unaltered from the commit f6288d6 and explicitly labelled. No empirical physical prediction or emergent geometry has been earned.

The next decisive **apparatus** experiment is a fixed mass × window-width matrix with independent clock positional record statistics and detector output, plus spatial-grid/time-step convergence for *each* point. The stronger **theory** hurdle remains a physically implemented source/detector coupling without privileged global access to the causal configuration, together with dynamical selection and a Lorentz-compatible propagation law. A healthy finite clock cannot, by itself, earn either of these.

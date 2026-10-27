# Field action, particles and coherent matching

The active equations and code units are unchanged from 7d9c93b. With three real canonical fields \(\chi=(\chi_1,\chi_2,\chi_\phi)\), momenta \(\pi=\partial_s\chi\), and \(r^2=\chi_1^2+\chi_2^2\),

\[
\mathcal H=\frac12\pi^2+\frac12(\nabla\chi)^2+
\frac{\lambda_X}{4}r^4+\frac{g_X^2}{2}\chi_\phi^2r^2+
\frac{\lambda_\phi}{4}\chi_\phi^4,
\qquad (\lambda_X,g_X^2,\lambda_\phi)=(0.001,0.01,0.03).
\]

The fixed radiation clock is \(\tau=20+a_m s+Cs^2/2\), \(a_m=8.869218738806\), \(C=1.9857962750634526\). The spatial derivative is the reviewed FD4 operator. Its nonnegative symbol is used consistently in the energy, tangent evolution and exports. The shared stepping interface differentiates the actual discrete Yoshida map; it does not evolve a different tangent approximation.

## What is exactly conserved

The periodic Hamiltonian conserves its energy and

\[
Q=\int d^3x\,(\chi_1\pi_2-\chi_2\pi_1).
\]

There is no independent exact symmetry conserving particle plus antiparticle number in this relativistic quartic theory. At any fixed positive reference \(\omega_k\),

\[
N_\omega=\sum_{k,a=1,2}\frac{|\pi_{a,k}|^2+\omega_k^2|\chi_{a,k}|^2}{2\omega_k},
\qquad N_{\pm,\omega}=\tfrac12(N_\omega\pm Q).
\]

Writing \(\dot\pi_{a,k}=-\omega_k^2\chi_{a,k}-I_{a,k}\), the exact identity is

\[
\dot N_\omega=-\sum_{k,a}\frac{\operatorname{Re}(\pi_{a,k}^*I_{a,k})}{\omega_k}.
\]

A time-dependent reference adds its explicit derivative. These reference changes are not by themselves a physical chemical rate. Likewise, observed near constancy of a frequency-weighted field statistic establishes approximate action transport over the measured interval; it is not a new Noether law or a bound on the much later binding epoch.

## The measured negative-frequency statistic

The complex classical field is \(\psi=\chi_1+i\chi_2\), so a lab-frame \(e^{+i\omega s}\) component carries positive current in this convention. For a Hann-tapered window of \(M\) uniformly spaced samples,

\[
P_k(\omega_n)=\frac{|\sum_j h_j\psi_k(s_j)e^{-i\omega_ns_j}|^2}{M\sum_j h_j^2},
\qquad \mathcal A_- = \sum_{k,\omega_n<0}|\omega_n|P_k(\omega_n).
\]

The coherent zero mode and connected nonzero modes are reported separately. In a free stationary theory this reproduces the usual occupation, before the stated quantum vacuum convention. In an interacting coherent medium it is a well-defined classical frequency action. A small spread under different positive reference masses makes it a useful diagnostic; that property alone does not specify a quantum particle basis or a Boltzmann distribution.

The spectral arrays retain \(k\le8\). The exact equal-time reference contribution above that band is reported separately, rather than silently discarded. Frequency resolution is set by the window length; zero padding would not improve it. Statistical spectra and the full matrix retarded response are different inputs. Failure of an equilibrium fluctuation-dissipation relation alone does not rule out nonequilibrium quasiparticles. A kinetic assignment still needs justified mode residues, mixing, widths and coherent separation where most of the population sits.

## Why the old infrared vacuum cannot just be extended

At the initial preparation slice, the six-real-field conformal Hessian has eigenvalues about \(-2.83238,-1.99965,15.6815,529.063,1320.08,3201.44\). Thus \(k^2+\lambda_a\) is negative for two directions below approximately 1.683 and 1.414. An inverted oscillator has no instantaneous ground state. This is a state-preparation obstruction, not a numerical failure or a prohibition on Gaussian states.

For the added modes, each stable direction retains \(\Omega_a=\sqrt{k^2+\lambda_a}\). In a negative-curvature direction the two controls supply \(\Omega_a=\sqrt{k^2+\sigma^2}\), with \(\sigma=1,3\). In the initial Hessian eigenbasis,

\[
C_{qq}=\tfrac12\operatorname{diag}(\Omega_a^{-1}),\qquad
C_{pp}=\tfrac12\operatorname{diag}(\Omega_a),\qquad C_{qp}=0.
\]

These are allowed pure canonical Gaussian reference states of the unconstrained scalar system. They are not ground states, a full gauged state or a selected inflationary covariance. The full scalar transfer and independent neutral transfer carry them to \(\tau=20\); only the X/phi marginal is subsequently evolved nonlinearly. The old seeded modes are left unchanged. No sampled charge or energy is retuned.

## What coherent backward rotation means

Once a free massive description applies, \(\psi=Ae^{+ims}+Be^{-ims}\) has

\[
n_+=m|A|^2,\qquad n_-=m|B|^2,\qquad Q/V=n_+-n_-.
\]

The backward-rotating expectation value therefore contributes to antiparticle number in that basis, including for a coherent quantum state. It is not an incoherent thermal gas. A separate postulate deciding whether to count it is unnecessary. The unresolved input is the dynamics that connects the interacting coherent field to that particle basis.

For a homogeneous central potential, let \(q=Q/V>0\),

\[
H_r=\frac{p_r^2}{2}+\frac{q^2}{2r^2}+\frac{m^2r^2}{2}+\frac{\lambda r^4}{4},\qquad
J_r=\frac1\pi\int_{r_-}^{r_+}dr\,\sqrt{2E-q^2/r^2-m^2r^2-\lambda r^4/2}.
\]

In a free isotropic harmonic oscillator, \(E/m=2J_r+q\), hence \(n_-=J_r\). A slow central deformation preserves \(J_r\) and charge. The computed matching control uses \(m(s)=30 f(s/T)\), \(\lambda(s)=\lambda_X[1-f(s/T)]\), \(f(u)=[1-\cos(\pi u)]/2\). It ends in a free theory and explicitly checks work as well as charge. A sudden deformation supplies a different amount of work at different radial phases and need not give the same final number. Neither deformation is claimed as the hidden sector's binding history.

For the reviewed starting orbit, the lab-frequency statistic and radial action are close but distinct. Their numerical values and both the coherent-charge and original full-state-net-charge denominators are provided in COHERENT_MATCHING.json. The approximately 1.74 statistic can be retained as a leading classical action budget; it cannot be promoted to a preparation-independent, exactly conserved out-particle population.

## Averaged export

Each export retains \(\langle C_{ij}(k)\rangle\), its FD4-\(k^2\)-weighted counterpart, \(\langle v_{0,i}v_{0,j}\rangle\) for the coherent mode, and the averaged full potential. Averaging the mean vector and squaring it would erase coherent kinetic/current contributions, so it is not used. Both uniform-conformal and uniform-cosmic-time averages are supplied with their windows and exact charge/energy reconstructions. The scalar/interaction energy convention is unchanged.

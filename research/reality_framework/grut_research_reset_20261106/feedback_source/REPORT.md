# Keep the fermion source paused and advance a testable cogenesis comparison

The refined one-mode force missed the registered feedback gate on
the transverse cutoff. Changing the cutoff from $6M$ to $8M$ shifted
the instantaneous force by **4.36% of peak**, against a 3% limit.
The same change moved the net work by only $6.30\times10^{-6}$
relative. Keep this source as a prescribed-background calculation.
The live holdout stays unopened.

Prioritize axion-like-particle (ALP) cogenesis as the competing
matter mechanism for observational work. A fixed electroweak
anomaly texture and a specified charge-to-particle transfer yield
turn the dark and baryon densities into a microwave power prediction.
A published haloscope already has sensitivity to this comparison.
The transfer yield and anomaly texture are the conditions to test.

[First protocol](PROTOCOL.md) · [Registered refinement](REFINEMENT.md) ·
[Force derivation](METHOD.md) · [Refined gate](results/refinement/prescribed_gate.json) ·
[Cogenesis calculation](results/cogenesis_screen.json)

## What the source test resolved

Restore the full coherent force, keep a constant initial-vacuum
energy subtraction, and fix the finite quadratic terms from the
retarded continuum response. The former pointwise subtraction
changes the instantaneous force by a conservative term whose
closed-pulse integral vanishes. At the refined grid it gives
a peak of $1.510M^4$, compared with $0.02436M^4$ for the
matched source. The finite kinetic matching is made for the
selected Fourier mode with an invariant continuum reference.

The first registered gate failed. Its fine matched force disagreed
with full retarded linear replacement by 53.25% of peak. The
second-order time split shifts the static susceptibility by
$0.01231M^4$ at that grid. Its leading error scales as
$\Delta t^2\Lambda^4$ for a common large momentum cutoff.
Correcting that linear artifact reduced the discrepancy to about
1.50%, while a separate coarse-momentum shift of about 42% remained.
[First gate](results/prescribed_gate.json) ·
[Time-split diagnosis](results/split_diagnosis.json)

The fourth-order control showed force error ratios of 16.51 and
16.06 on successive timestep halvings. Its folded integral
agreed with the full momentum integral to $2.8\times10^{-15}$
relative. The retarded kernel qualified against the parameter
integral and the weak spectral work before the source comparisons.

| Refined comparison | Force difference as % of fine peak | 3% gate |
| --- | ---: | --- |
| Halved timestep | 0.29 | Pass |
| Denser Bloch and radial quadrature | 1.34 | Pass |
| Increased longitudinal band count | 0.58 | Pass |
| Transverse cutoff $6M$ to $8M$ | **4.36** | **Miss** |
| Base to combined fine grid | 1.34 | Pass |
| Fine open to tapered boundary | 0.28 | Pass |
| Fine local matching to full linear replacement | 2.72 | Pass |

Full linear replacement leaves about 4.0% of peak in the
radial-cutoff comparison. This residual lies in the nonlinear
force in this decomposition. Fine net work is $0.134863972M^4$
and free outgoing energy is $0.134863770M^4$, with a relative
difference of $1.50\times10^{-6}$. The remaining nonlinear momentum tail requires its own bound. Reopening feedback would require a
controlled transverse tail estimate or an independently converged
force on a larger cutoff ladder. The test uses flat space, one
scalar mode and one Majorana species, with a specified vacuum.
Cosmological feedback would also require the full spatial stress
and matching to the dynamical background.

## The observational decision

H55 still requires a physical spatial source, a reciprocal
stress and a controlled transport transition before it predicts
an independent observable. M1 has a conditional multi-target
recoil shape, but its archived portal-heating budget gives
about 6–22 events at the assumed exposure and the recoil node
is much weaker than the peak. ALP cogenesis offers a direct
comparison using two cosmic densities as calibration and a
haloscope coupling as an independent measured quantity.
[Prior observable screen](../strong_drive/OBSERVABLES.md) ·
[Co, Hall and Harigaya](https://doi.org/10.1007/JHEP01(2021)172) ·
[ORGAN Phase 1b](https://doi.org/10.1103/PhysRevLett.132.031601)

For the comparator, use the leading-order electroweak EFT with
only the $SU(2)_L$ anomaly coefficient $c_W=-1$ in the authors'
convention, with $c_Y=c_g=0$ and other ALP couplings set to zero
at the electroweak matching scale. The positive rotation then
biases electroweak sphalerons with $c_B=12/79$ and the photon
coefficient is $c_\gamma=c_W+c_Y=-1$ at leading order.
Assume $f_a(T_{\rm EW})=f_a$, rotation through sphaleron freeze-out,
and conserved angular yield $Y_\theta$ between freeze-out and
conversion. In the kinetic-misalignment branch the transfer
$C=Y_a/Y_\theta$ gives

$$
Y_B=c_B Y_\theta\frac{T_{\rm EW}^2}{f_a^2},\qquad
\frac{\rho_a}{s_{\rm ent}}=C\,m_aY_\theta,\qquad
|g_{a\gamma\gamma}|=\frac{\alpha|c_\gamma|}{2\pi f_a}.
$$

The unknown initial charge cancels from the abundance ratio.
Fixing $m_a$, the texture and $C$ then fixes $f_a$ and the photon
coupling using the two measured densities. Allowing an independent
hypercharge anomaly, for instance $c_Y=-c_W$, can cancel the
leading photon signal while keeping the sphaleron coefficient.
A fixed texture is necessary to predict a photon signal.

The haloscope power contains $g^2\rho_{a,\rm loc}/m_a$. Let
$f_{\rm cos}$ and $f_{\rm loc}$ be the ALP fractions of cosmic
and local dark matter and write $d=\rho_{\rm DM}/s_{\rm ent}$.
In this model,

$$
\frac{g^2\rho_{a,\rm loc}}{m_a}
=\left(\frac{\alpha c_\gamma}{2\pi}\right)^2
 \frac{C Y_B\,\rho_{\rm DM,loc}}{c_B T_{\rm EW}^2 d}\,
 \frac{f_{\rm loc}}{f_{\rm cos}}.
$$

With fixed transfer and a population that traces the local dark
matter, reducing the cosmic fraction increases the coupling
in the abundance matching and leaves this power prefactor
unchanged. Detector response still depends on frequency.
Fragmentation and halo substructure can change $C$ and the
local-to-cosmic fraction ratio.

As a retrospective example choose $m_a=110\,\mu\mathrm{eV}$
inside the published ORGAN window. Use $C=2$, $T_{\rm EW}=130$ GeV,
$Y_B=8.7\times10^{-11}$, $d=0.44$ eV and all dark matter in the
ALP. The calculation gives $Y_\theta=2000$, $f_a=2.43\times10^8$
GeV, and $|g_{a\gamma\gamma}|=4.78\times10^{-12}\,\mathrm{GeV}^{-1}$
at 26.60 GHz. At freeze-out the estimated rotation rate is
$3.17\times10^7 m_a$ and its energy is $3.58\times10^{-5}$
of radiation, consistent with the assumed kinetic and
radiation-dominated limits at that temperature.

ORGAN reports sensitivity of $4\times10^{-12}\,\mathrm{GeV}^{-1}$
over 107.42–111.93 $\mu\mathrm{eV}$ at 95% confidence for its
standard local halo density of $0.45\,\mathrm{GeV/cm^3}$ and
all of it in axions, with a loss of coverage near 26 GHz.
The example lies away from that notch and gives 1.43 times
the power at the quoted coupling sensitivity. This is a
screen against the published result, with the example mass
chosen after the window was known. At the quoted sensitivity the threshold in this example
is approximately $C f_{\rm loc}/f_{\rm cos}=1.40$ with the
texture and other inputs fixed.
[ALP cogenesis, Eqs. 2.8, 2.9, 2.14 and 3.10](https://arxiv.org/abs/2006.04809) ·
[ORGAN, Eq. 1 and Fig. 5](https://doi.org/10.1103/PhysRevLett.132.031601)

The next discriminating calculation is to fix the anomaly
texture and preparation and determine $C$ through trapping
and fragmentation, independently of the haloscope data.
Use the resulting coupling and its uncertainty in an unused
frequency window. The rotating condensate can fragment, so
using $C=2$ from homogeneous trapping as an exact input would
prematurely fix the key uncertainty. The relation and existing
test belong to published cogenesis physics. They provide a
comparator for GRUT's search for a constraining mechanism.
[Kinetic fragmentation](https://arxiv.org/abs/2206.14259)

## Replay

Run from this directory with Python, NumPy and SciPy. The
prescribed-gate commands below return a nonzero status on the
archived misses. The data under [results](results/) include
complete sampled source histories and the gate metrics.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python source.py qualify
for case in base time spacing longitudinal transverse fine tapered; do
    python source.py "$case"
done
python gate.py prescribed
python diagnose.py
python higher_order.py qualify
for case in base time spacing longitudinal transverse fine tapered; do
    python higher_order.py "$case"
done
python higher_order.py gate
python cogenesis_screen.py
```

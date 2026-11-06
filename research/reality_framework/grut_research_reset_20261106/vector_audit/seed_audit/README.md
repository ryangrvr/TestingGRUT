# Radiation seed audit for the frozen 1 meV vector benchmark

**The frozen conversion fails as an ordinary local radiation gas.**
At the end of either inflation blend the assigned radiation energy is
far below the energy of a single horizon-scale quantum, even in the
initial volume that can causally feed a point at the nominal transition.
The logistic closure therefore remains a prescribed effective history.
Its physical seed, energy transfer, and perturbative preparation remain
open. The replicated common-spectrum abundance enhancement of 4.95
belongs to the specified expansion histories. This audit holds the mass,
inflation parameters, transfer result, and fitted \(r_{0.05}=0.002765\)
fixed. [Frozen replication](../replication/README.md).

## The required initial inventory

Let \(N=0\) mark the reference end of inflation. On the stiff solution
\(w_s=3/5\), \(N_i=\delta\) is the end of the inflation blend and
\(N_R=25.0037\) is the nominal half-radiation point. The local law
\[
Q=(\sigma-0.8)H\,\rho_s\rho_r/\rho
\]
gives \(df/dN=\sigma f(1-f)\) for \(f=\rho_r/\rho\). Thus
\(\ell_i=\log[f_i/(1-f_i)]=\sigma(N_i-N_R)\).
The zero-seed solution is exactly \(f(N)=0\) by uniqueness. Starting
from a vacuum bath requires an initial-state preparation or a source
that operates at \(f=0\). For a local transfer \(S\) from the stiff
component into radiation, \(df/dN=\sigma f(1-f)+S/(H\rho)\).
In the small-\(f\) regime, the weighted early source
\(\int dN\,e^{-\sigma(N-N_i)} S/(H\rho)\) adds directly to the
effective initial seed. The origin and timing of a seed source matter.

| Inventory at \(N_i\) | Narrow, \(\sigma=24\) | Wide, \(\sigma=12\) |
|---|---:|---:|
| \(f_i\) | \(1.02\times10^{-260}\) | \(2.98\times10^{-130}\) |
| \(H_i\) in GeV | \(5.61\times10^{12}\) | \(4.52\times10^{12}\) |
| \(\rho_{r,i}=3\bar M^2 H_i^2 f_i\) in \(\mathrm{GeV}^4\) | \(5.72\times10^{-198}\) | \(1.08\times10^{-67}\) |
| Hubble radius in metres | \(3.52\times10^{-29}\) | \(4.37\times10^{-29}\) |
| Energy in a Hubble sphere in GeV | \(1.36\times10^{-235}\) | \(4.91\times10^{-105}\) |
| Upper count of quanta with energy \(\ge H_i\) in that sphere | \(2.42\times10^{-248}\) | \(1.09\times10^{-117}\) |
| Upper count in the initial causal sphere of a point at \(N_R\) | \(2.78\times10^{-203}\) | \(8.55\times10^{-73}\) |
| Upper count in a pivot sphere of radius \(a_i/k_p\) | \(2.46\times10^{-170}\) | \(7.56\times10^{-40}\) |
| Formal ideal-bath \(T_i/H_i\) | \(3.58\times10^{-63}\) | \(1.65\times10^{-30}\) |

For a kinetic radiation gas, a quantum localized on the Hubble scale
has physical energy at least of order \(H_i\). Counting at that lower
energy gives \(E_i/H_i=4\pi(\bar M/H_i)^2 f_i\) as the number bound.
For the larger causal sphere use the initial radius
\[
 R_i\le H_i^{-1}
    \frac{\exp[(2.4-1)(N_R-N_i)]-1}{2.4-1}.
\]
Pure stiff evolution gives this upper bound because adding radiation
slows the drop in \(H\). This permits propagation at the speed of light
through the nominal wait. A nonnegative particle count with this mean
bounds the probability of even one such seed quantum by the same
number. One horizon-energy quantum per Hubble sphere would instead
give \(f_i=4.22\times10^{-13}\) or \(2.74\times10^{-13}\). Under the
homogeneous law this would give the midpoint at \(N=1.25\) or \(2.56\).

The temperature row uses the archive's ideal \(g_*=106.75\) solely
as a density conversion. The microscopic wavelengths associated
with that formal temperature far exceed the Hubble radius. A vacuum
or coherent long-wavelength state would need its own stress tensor
and coupling. The count argument applies to a causal gas of on-shell
quanta.

## Which radiation is produced

The frozen quadratic calculation already produces relativistic
longitudinal Proca modes and both tensor helicities. At \(N=4\), the
resolved \(8\le K=k/k_{\rm end}\le16\) band alone gives longitudinal
energy fractions \(3.11\times10^{-13}\) and \(3.96\times10^{-13}\).
The two tensor helicities give twice those values. The archive keeps
these as spectator stress. This audit uses the frozen hard occupations,
with the time-resolved hard-pressure subtraction in the archived audit.
A law that included every relativistic
species in \(f\) would receive a much larger seed. A visible-bath law
requires specifying which species participate. [Vector production
mechanism](https://arxiv.org/abs/1504.02102).

At the free-field level, the massless Maxwell action in four-dimensional
FLRW reduces in conformal coordinates to the flat-space action.
A conformal Maxwell vacuum therefore stays empty under homogeneous
expansion at that order. The transfer fit supplies an ideal thermal
bath normalization. Curvature couplings, a visible-field initial
state, and the interaction responsible for populating that bath are
additional inputs. The published Higgs-production rates from an
oscillating inflaton condensate use a particular background and
couplings. [Gravitational portals](https://arxiv.org/abs/2112.15214).

A concrete gravitational source appears in a completion with one
initially empty, minimally gravitating Maxwell field. Deep inside
the horizon the relativistic longitudinal vector has the
Stückelberg scalar limit. Two of these quanta can create a Maxwell
pair through Einstein gravity even with zero visible occupancy.
For that one massless gauge field the relativistic cross section
is \(s/(480\pi\bar M^4)\), including the identical final-state
factor. The retained momenta have \(p/H>2160\) and
\(p/\bar M<7.8\times10^{-7}\) at \(N=4\).
[Gravitational cross sections, Eqs. 3.7 and 3.12](https://arxiv.org/abs/1708.05138).

For a phase-averaged isotropic kinetic estimate, restrict incident
angle cosines to \([-1/2,1/2]\). This avoids the exactly opposite
squeezed partner directions. For massless quanta
\(s=2p_1p_2(1-\mu)\) and the relative flux factor is \(1-\mu\).
The angular integral is \(13/24\). With one longitudinal degree
of freedom, write \(\rho_L=\int p\,dn_L\) and \(J_2=\int p^2\,dn_L\)
in the resolved band. The rate of events per physical volume is
\((13/24)\rho_L^2/(480\pi\bar M^4)\) and the injected power is
\(2(13/24)\rho_L J_2/(480\pi\bar M^4)\).
Integrate the source and redshift over \(4\le N\le4.1\) before
applying any \(Q\) amplification. This uses a kinetic approximation
to the produced state. The quadrature changes the injection by
less than 0.04% across the archived meshes.

| One Maxwell field example | Narrow | Wide |
|---|---:|---:|
| Direct injected energy fraction at \(N=4.1\) | \(2.02\times10^{-45}\) | \(3.24\times10^{-45}\) |
| Excess over the frozen seed there, orders of magnitude | 173 | 64 |
| Expected events in one initial comoving Hubble sphere over the interval | \(1.15\times10^{-28}\) | \(1.87\times10^{-28}\) |
| Midpoint if injected energy immediately behaved as homogeneous \(Q\) radiation | 8.39 | 12.64 |

This example produces a mean that overwhelms the frozen visible
seed if that energy feeds the prescribed fluid. The small number
of events per Hubble sphere prevents a homogeneous-fluid
interpretation at production. An actual conversion time would
depend on which particles feed the bath, their spatial distribution,
and the dynamics of amplification and thermalization. The estimated
Maxwell rate is a channel in the stated completion with the stated
kinetic approximation.

## Independent seed perturbations and the scalar calibration

Start on a uniform-total-density slice and vary \(\ell_i\)
independently of the inflaton mode. With \(n=N-N_i\), the exact
homogeneous family gives
\[
 \log(\rho/\rho_i)=-4.8n+
 \frac{0.8}{\sigma}
 \left[\log(1+e^{\ell_i+\sigma n})-\log(1+e^{\ell_i})\right].
\]
At late fixed total density,
\(dN/d\ell_i=(0.2/\sigma)(1-f_i)\).
Direct finite shifts recover this derivative for both histories.
In the long-wavelength limit let \(b\) be the root power per
logarithmic wavenumber of the initial independent logit
perturbation. The added curvature power is
\((dN/d\ell_i)^2 b^2\) for zero cross correlation with the
inflaton mode. [Adiabatic and entropy transfer](https://arxiv.org/abs/astro-ph/0211602).

| \(b\) at the frozen pivot | Added power divided by \(A_s\), narrow | Added power divided by \(A_s\), wide |
|---|---:|---:|
| \(10^{-4}\) | \(3.31\times10^{-4}\) | \(1.32\times10^{-3}\) |
| \(10^{-3}\) | 0.0331 | 0.132 |
| \(10^{-2}\) | 3.31 | 13.2 |

A 10% budget relative to the frozen \(A_s=2.1\times10^{-9}\)
requires \(b<1.74\times10^{-3}\) or \(b<8.69\times10^{-4}\).
These are conditional sensitivity thresholds. A seed state would
have to supply its own perturbation spectrum and correlations.
The common time-shift solution alone leaves this independent mode
free. The reported \(r_{0.05}=0.002765\) retains the original
adiabatic scalar calibration.

## Reproduction and standalone validation

Use the same NumPy and SciPy versions as the replication packet.
Run this audit from the repository with the adjacent frozen
replication input files and mode spectra.

```sh
python audit.py
```

[Inputs and protocol](PROTOCOL.md). [Computed values](results.json).

For use outside the repository, set `GRUT_REPLICATION_PACKET` to the
directory containing the extracted ZIP's `frozen_inputs.json` before
running the same command.

The distributed ZIP already contains `python replicate.py "$history"`.
The incorrect string appears in the pasted command. The archived
README needed no file edit. In a fresh Python 3.12.14 environment,
the command block extracted verbatim installed the pinned NumPy
2.3.5 and SciPy 1.17.0 dependencies. From the extracted ZIP, all
eight runs and both comparison ledgers reproduced byte for byte.
[Validation record](standalone_validation.json).

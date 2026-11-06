# X11 · Finite H55 evolution and matching equations

## Shared-metric production

Start with the archived H55 fields and velocities at $tH_i=-0.02$.
Project the same 64³ Gaussian realization onto each centered Fourier
cube. Evolve three scalar fields, one common scale factor, and the four
sterile species at each cell. Each cell carries both helicities and the
radial shell weights of X5. The initial Bloch vectors use the X5 dressed
vacuum, with the local scalar gradient included in the initial classical
acceleration. This is a specified state in the cell-local model.

The cell action uses X5's inertia $A_j(a)=1-Z_j(a)$ and its matched
gravitational constraint. The same $A_j$ multiplies the spatial gradient
energy, as in X10. At fixed scale, for spatial step $d$ the added comoving
energy is
$$
E_G = \frac{a}{2d^2}\left\langle
\sum_{j,r} A_j(a)[\phi_j(\mathbf x+d\mathbf e_r)-\phi_j(\mathbf x)]^2
\right\rangle .
$$
The sum uses the three forward links. The force is $a A_j\nabla_d^2\phi_j$
and the associated pressure is $-(A_j+\partial_{\ln a}A_j)G_j/(3A_j)$.
The scalar site-to-site kick brackets the X5 midpoint exchange. At
fixed volume the X5 rotation supplies the same canonical impulse to
the scalar that changes the sterile energy. The separately summed
comoving quantum work uses $-\bar v_j\Delta\Pi_{j,N}$ and the
quantum volume work uses the quantum stress. The total stress and
Friedmann constraint use volume averages. The recorded phase
operator deviations use the actual complex phases of the cell fields.
[X5 local action](../x5_coupled/MODEL.md) · [X10 gradient matching](../x10_handoff/MODEL.md)

The fermion continua lack links between sites. The ultraviolet
subtraction is inherited from the local homogeneous theory. Spatial
fermion propagators, bath covariances, finite-density effects and
covariant matching of the gradient terms are physical extensions.
Real-time fermion ensemble methods give one route to propagating
a field with quantum backreaction on a lattice.
[Borsányi and Hindmarsh](https://doi.org/10.1103/PhysRevD.79.065010)

## Conditional kinetic and release model

The projection onto the local free mass axis gives $0\le f_{hj}\le1$
and preserves the free mass number and energy. Averaging those
occupations across the sites is linear in both quantities. The
kinetic run keeps all native shell nodes and weights and appends
empty shells at higher comoving momentum. The two scalar parents
are assigned harmonic reservoir energies in proportion to their
positive harmonic diagnostics at the boundary, normalized to their
represented field energy. The scalar derivative subtraction already
enters that energy. The spectator density and its dark radiation
are carried separately. The pressure discontinuity and the omitted
coherent scalar work are measured before applying this approximation.
[X4 kinetic ledger](../x4_closure/REPORT.md)

A second-order expansion step uses $\dot H=-(\rho+P)/(2M_P^2)$
after the fixed-volume transfers. The release solver uses a single
cubic interpolation of $\ln a$ with the supplied Hubble slopes.
Its first and second derivatives supply $H$ and $\dot H$. The
background trace input is the cycle-averaged production stress and
the pressure trace of the conditional kinetic model. The spectator
mass uses the evolved $g_s^2\langle\chi^2\rangle$ before the
handoff and the inherited harmonic form afterward. The conditional
release uses the M2-X holder product $c_0=1.55$, zero holder offset
and hidden coupling $e=10^{-6}$. The product is implemented in the
holder coupling with the spectator amplitude fixed at the archived
H55 value. The H55 scale stays fixed.

The visible entropy is inferred from the final visible radiation,
with the stated SM degrees of freedom and no later visible injection.
The produced dimensionless comoving charge is $N_X$. Using the same
normalization for scale factor and fields gives
$$
Y_Q = \frac{|N_X|}{a^3 s_{\rm vis}/H_i^3},\qquad
\Omega_X h^2 = \frac{s_0}{\rho_{\rm crit}/h^2}\,M\,m_X Y_Q .
$$
The inherited population map receives this calculated $Y_Q$.
$M$ is its final mass per unit charge in units of $m_X$.
The Q-matter and free-X prescriptions carry their original
approximation and extrapolation flags. The comparison target is
$\Omega_c h^2=0.120$ under base $\Lambda$CDM.
[Fixed-yield population map](../e2e_relic/code/relic.py) ·
[Planck cosmological parameters](https://arxiv.org/abs/1807.06209) ·
[PDG cosmological constants](https://pdgweb.lbl.gov/2024/download/db2024.pdf)

## Axial CP matching postulate

In the specified inflaton interaction basis, postulate that
the sterile Yukawa columns share one real flavor convention.
The physical Majorana masses at the vacuum are
$\eta_j M_j\exp(i\theta_j)$, with $\eta_j=\pm1$
a discrete parity and $\theta_j=\phi_{j,\rm vac}/f_j$.
The literal law has all $\eta_j=+1$. The extension uses the
inherited parity $\eta_2=-1$. In the positive mass basis,
$Y_m=Y_R\operatorname{diag}(d_j\exp(-i\theta_j/2))$
with $d_j^2=\eta_j$. The decay prediction is
$$
\frac{\operatorname{Im}[(Y_m^\dagger Y_m)_{ij}^2]}
     {|(Y_m^\dagger Y_m)_{ij}|^2}
  =\eta_i\eta_j\sin(\theta_i-\theta_j)
$$
for nonzero overlaps. In H55, the shared waterfall mass phase
sets all three normalized invariants within N2–N4 to zero.
The other three values follow from the one relative axial angle.
The unflavoured asymmetries use the inherited vertex and
self-energy loop functions.
[Covi, Roulet and Vissani](https://arxiv.org/abs/hep-ph/9605319)

Let $u_j=v_{\rm EW}Y_{R,j}/\sqrt{M_j}$ so that
$u_j^T u_j=\widetilde m_j$ fixes the total decay width.
The light mass matrix at tree matching is
$m_\nu=-\sum_j \eta_j e^{-i\theta_j}u_j u_j^T$.
Its Takagi mass sum is its nuclear norm and obeys
$$
\sum_a m_{\nu,a}=\|m_\nu\|_*\ge |\mathrm{Tr}\,m_\nu|
 =\left|\sum_j\eta_j e^{-i\theta_j}\widetilde m_j\right| .
$$
The bound holds for all real column shapes at the supplied
widths. Oscillation and loop comparisons for the parity
extension use the fixed benchmark texture. The late baryon
comparison carries the corresponding widths into the
kinetic history and uses the inherited unflavoured network
with its thermal equilibrium subtraction. It is a conditional
test of the matching postulate on those inputs.

The postulate adds information beyond closed-system dynamics.
For fixed positive masses, multiplying Yukawa columns by phases
leaves every total sterile width invariant while changing these
CP-odd invariants. The collisionless source is independent of the
Yukawa matrix. Production and total widths therefore cannot
select this matching rule on their own.

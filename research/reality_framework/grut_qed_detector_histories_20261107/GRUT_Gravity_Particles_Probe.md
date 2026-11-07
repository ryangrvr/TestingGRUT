# Can GRUT get gravity and particles?

A probe of the GRUT integrated candidate at [TestingGRUT commit 8e90c36](https://github.com/ryangrvr/TestingGRUT/tree/8e90c36ecf8f0986fbccb1097b327f37526959c4), 6 November 2026.

This is an outside probe, not part of the GRUT record. Nothing was committed to the repository. The scripts `particles.py` and `gravity.py` reproduce every number below in under a minute (NumPy and SciPy only).

## Short answer

- **Particles: partly.** The photon bridge's recipe (a local lattice Hamiltonian on a spatial torus with continuous time) can host photons, and with standard lattice techniques it can also host vector-like charged matter such as QED electrons and QCD quarks. It cannot host the Standard Model's chiral weak sector without new structure. Probe P1–P2 shows why: every lattice of this type produces left- and right-handed particles in equal numbers.
- **Gravity: only as adopted linearized content.** A spin-2 field can be placed on the same lattice. It gives Newton's law, two graviton polarizations and the correct light bending, and the observed light bending rules out the obvious alternatives. But the coupling to matter has to be supplied, which matches the repository's own September finding. Full nonlinear gravity with dynamical geometry is not reachable from a fixed lattice without a further principle.
- **New constraint on the photon bridge.** The lattice's own photon dispersion, compared with LHAASO gamma-ray burst timing, limits the regulator spacing to about \(a \lesssim 10^{-27}\) m. That gives the photon bridge's "physical spacing must be identified" obligation a numerical ceiling.

## What the repository already says

| Source in the repo | Status of gravity and matter |
|---|---|
| `GRUT_Integrated_Theory.md`, Section 7 | Semiclassical GR plus the Standard Model is the adopted target, eqs. (12)–(13). Matching to it, eq. (14), is an explicit unmatched bridge. "GR and the Standard Model enter as adopted content." |
| `PHOTON_MATCHING.md` | Free Maxwell photon on a fixed flat torus. Charged matter and metric response are omitted. |
| `GR1_GRAVITY_VERDICT_01.md` (25 Sep) | Gravity's spectral exponent is "selected by structure in class", but the minimal-stress coupling and a gapless retained sector are **load-bearing imports**. "No Einstein equations … were used or implied." |
| `books/BOOK_VII_MATTER_SM_INTERFACE.md` | Self-described "absence map": most of the matter sector has no GRUT account on the record. |
| `grut_research_reset_20261106/AUDIT.md` | Fermion production work (H55/X11) takes the Majorana particle dictionary, the FRW metric and visible-sector laws as inputs. Its handoff fails pressure and work criteria. |

So the record is consistent: neither gravity nor Standard Model particles are derived. The probes below ask what the current construction method could support next, and where it gets stuck.

## Particles

### P1. A naive lattice fermion comes in eight copies

A single right-handed Weyl fermion discretized like the photon, \(H(k)=\sum_i \sigma_i \sin k_i\), has zero-energy nodes at all eight corners \(k\in\{0,\pi\}^3\) of the Brillouin zone.

| Species | Right-handed | Left-handed | Net chirality |
|---:|---:|---:|---:|
| 8 | 4 | 4 | 0 |

### P2. Random local changes don't fix it

Adding random Hermitian hopping terms up to two sites away, and locating every node with a root search and its Jacobian sign, gives:

| Deformation scale | Nodes found | Right | Left | Net |
|---:|---:|---:|---:|---:|
| 0.03 | 8 | 4 | 4 | 0 |
| 0.08 | 8 | 4 | 4 | 0 |
| 0.15 | 8 | 4 | 4 | 0 |
| 0.25 | 8 | 4 | 4 | 0 |

An independent Berry-flux sweep (the Chern number of each \(k_z\) slice) changes exactly at those nodes. This is the [Nielsen–Ninomiya theorem](https://en.wikipedia.org/wiki/Nielsen%E2%80%93Ninomiya_theorem) at work: a local, Hermitian, translation-invariant lattice with conserved charges has equal numbers of left- and right-handed fermions. In the electroweak sector, left- and right-handed fermions carry different charges, so this blocks a direct lattice version of the Standard Model.

### P3. Wilson fermions: one light particle, but no chiral symmetry

Adding a Wilson term with \(r=1\) and bare mass \(m=0.05\) (lattice units) gives one light Dirac particle. The seven doublers get masses of 2.05, 4.05 or 6.05 in units of \(\hbar c/a\). The light particle's dispersion matches \(E^2=k^2+m^2\) to within 0.7% at \(ka=0.02\), and the agreement is nearly direction-independent. The cost is that the Wilson term breaks chiral symmetry explicitly (commutator norm 8.2 at a doubler).

That is enough for vector-like theories, and it is the standard route in lattice QCD. It is not enough for the chiral SU(2)×U(1) weak sector. Known escapes (Ginsparg–Wilson/overlap fermions, domain walls in an extra dimension, mirror fermions) each add structure. A [lattice gauge theory review](https://ar5iv.labs.arxiv.org/html/hep-ph/9911400) notes that the no-go theorem "has so far precluded all attempts to achieve a realistic lattice formulation of the electroweak sector."

**What this means for GRUT:** this is a second selection problem of the same kind the photon bridge found. The access table and lattice rules do not say which escape route nature uses, and a choice has to be supplied.

## Gravity

### G1. The lattice gives Newton's law

The static exchange between two masses on the cubic lattice is the lattice Laplacian Green function, computed with the Bessel-integral representation. As a check, \(G(0)=0.252731\) matches Watson's classical value.

| \(r/a\) | Axis: \(4\pi r G\) | \(r/a\) | Diagonal: \(4\pi r G\) |
|---:|---:|---:|---:|
| 1 | 1.0815 | 1.73 | 0.9485 |
| 4 | 1.0198 | 6.93 | 0.9965 |
| 16 | 1.0010 | 27.7 | 0.9998 |

It converges to \(1/r\), and the anisotropy dies away within a few dozen spacings.

### G2. Spin-2 on the lattice has two positive-norm polarizations

Using the linearized Einstein (Fierz–Pauli) operator with lattice momenta \(\hat k_i=2\sin(k_i a/2)/a\):

- Linearized gauge invariance is exact, with residuals of \(10^{-15}\) or smaller at axial, off-axis and near-zone-edge momenta.
- The null space is 4-dimensional off shell (gauge only) and 6-dimensional on shell (4 gauge plus 2 physical).
- At the pole, the residue on the 6-dimensional space of conserved sources has **2 positive and 0 negative** eigenvalues for massless spin-2. For massive spin-2 it has 5 positive.

### G3. Calibrate on Newton, test on light bending

This mirrors the photon bridge's calibration-then-holdout logic. Each candidate's coupling is tuned so it reproduces Newton's law for slow masses, and it is then asked to predict how much light bends.

| Candidate force (Newton-calibrated) | Light bending vs GR | Implied PPN \(\gamma\) | Distance from Cassini |
|---|---:|---:|---:|
| Massless spin-2 coupled to full \(T_{\mu\nu}\) | 1.00 | 1.0 | 0.9σ |
| Massive spin-2, \(m\to0\) | 0.75 | 0.5 | ~21,700σ |
| Scalar coupled to \(T^\mu{}_\mu\) | 0 | −1 | ~87,000σ |

The comparison uses [Cassini's \(\gamma-1=(2.1\pm2.3)\times10^{-5}\)](https://pubmed.ncbi.nlm.nih.gov/14508481/). The 3/4 factor is the known [van Dam–Veltman–Zakharov discontinuity](https://arxiv.org/html/hep-th/0011138), and the scalar fails because light's stress tensor is traceless. These are textbook results reproduced in the lattice setting, not new physics. Their value is that the observed light bending selects among candidates that all reproduce Newton.

### Why this is still not "getting gravity"

1. **The coupling is put in by hand.** Universal coupling \(h_{\mu\nu}T^{\mu\nu}\) is what makes G3 come out right. That is exactly the import the repository's GR1 verdict found to be load-bearing ("minimal-stress: LOAD-BEARING"). Nothing in the access postulates forces it yet.
2. **Nonlinear gravity needs diffeomorphism invariance.** G2 has only the linearized gauge symmetry on a fixed flat lattice. Dynamical geometry, which is in the integrated paper's target (13), is not reached.
3. **Emergent gravitons face a no-go.** The [Weinberg–Witten theorem](https://en.wikipedia.org/wiki/Weinberg%E2%80%93Witten_theorem) forbids massless spin-2 particles carrying a Lorentz-covariant conserved stress tensor. A lattice is not Lorentz invariant at the cutoff, so a loophole exists. Using it, though, means showing that Lorentz invariance re-emerges for every species, to observational precision.

## G4. How fine the lattice must be

The plaquette-only photon of `PHOTON_MATCHING.md` has \(\omega=(2c/a)\sqrt{\sum_i\sin^2(k_ia/2)}\), so high-energy photons travel slightly slower:

| Direction | \(\Delta v/c\) | Spacing ceiling \(a_{\max}\) |
|---|---:|---:|
| Axis | \(-(ka)^2/8\) | \(9.3\times10^{-28}\) m |
| Face diagonal | \(-(ka)^2/16\) | \(1.3\times10^{-27}\) m |
| Body diagonal | \(-(ka)^2/24\) | \(1.6\times10^{-27}\) m |

The ceiling uses the [LHAASO GRB 221009A limit](https://ui.adsabs.harvard.edu/abs/2024PhRvL.133g1501C/abstract) \(E_{\rm QG,2}>6\times10^{-8}E_{\rm Pl}\approx7.3\times10^{11}\) GeV, in the convention \(\Delta v/c=-\tfrac32(E/E_{\rm QG,2})^2\). Published limits assume an isotropic effect, so treat these ceilings as order-of-magnitude.

So a GRUT-style lattice is allowed, but only if its spacing is below about \(10^{-27}\) m. That is still about \(10^{8}\) Planck lengths, so the Planck scale itself is not excluded. At that spacing, Wilson doublers would sit above about \(10^{11}\) GeV, far beyond collider reach.

## Recommended next steps for the program

1. **Add a vector-like charged fermion (lattice QED) to the photon bridge.** It is the smallest step that puts interacting particles into GRUT's detector-history framework with known controls. It would test the matching condition (14) beyond free fields.
2. **State the chiral obligation explicitly**, as the photon review did for the plaquette-only choice. Which escape from Nielsen–Ninomiya does GRUT adopt, and can access structure select it?
3. **Treat universal spin-2 coupling as the gravity target to derive.** G3 shows it is what observation demands, and GR1 shows it is the remaining import. Deriving it from access or closure would be a substantive result. Failing that, it should be stated as a permanent postulate.
4. **Carry the spacing ceiling into the photon bridge.** It turns "identify a physical spacing" into the constraint \(a\lesssim10^{-27}\) m, together with the requirement that photons, fermions and any graviton share an emergent light cone.

## Follow-up

All four next steps were followed through in `GRUT_Play_It_Out.md`. The detector-history checks pass with interacting QED, the preregistered spectrum checks miss at the current solver size, and a shared light cone tightens the spacing ceiling to about \(3\times10^{-32}\) m.

## Scope

- These are free and linearized calculations on a fixed flat lattice, written independently of the GRUT codebase.
- They reproduce known theorems (Nielsen–Ninomiya, vDVZ, Weinberg–Witten) and standard lattice results inside GRUT's construction style. The new part is applying them to GRUT's specific obligations and turning the spacing question into a number.
- None of this has been reviewed by the GRUT author or independent experts.

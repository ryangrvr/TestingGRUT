# A proposed empirical bridge, with a calculable recoil shape

22 October 2026 · parent 7abab86 · **proposal, not adopted Candidate 1**

Claude's central review is right: the mass/depletion identity is useful
internal physics, but the frozen gravity-only action has no demonstrated
laboratory access. The new calculation supplies a concrete possible bridge
and identifies its costs. It does not establish an allowed parameter region,
change the frozen action or explain dark energy.

## The useful construction

Research a separately versioned Higgs portal, adding to the potential

\[
\Delta U=\frac{\kappa}{2}(\phi^2-v_\phi^2)
\left(H^\dagger H-\frac{v_H^2}{2}\right).
\]

At leading portal order it turns the existing depletion identity into
an ordinary-matter source,

\[
g_{Qh}(0)=-\frac{\kappa v_H E_Q}{2\lambda_\phi v_\phi^2},
\qquad
\sigma_{QN}(0)=\frac{\mu_{QN}^2}{\pi}
\left[\frac{\kappa f_Nm_NE_Q}{2\lambda_\phi v_\phi^2m_h^2}\right]^2.
\]

The same solved field profile fixes a finite-momentum form factor F_Q(q).
It is not an extra freely fitted radius or scattering coefficient. The
form factor can change a visible recoil spectrum, unlike the previous
hidden D readout. See [MATHEMATICS.md](MATHEMATICS.md) for both derivations,
normalization, operator definitions and limits.

The bridge is conventional Higgs-portal physics, credited alongside related
soliton-detection work in [SOURCES.md](SOURCES.md). It is not a new interaction
derived from GRUT, nor is any scalar field identified with the SM Higgs.

## One qualification to the review

The low-energy/background-derivative credit belongs in the record. But
dE/dv_phi=E/v_phi does **not** by itself equal the physical radial field's
static source in this nonlinear model. Defining s=phi-v_phi, the static
equation gives

\[
J_s(0)=E/v_\phi+\lambda_\phi v_\phi\int s^2.
\]

The ratios J_s(0)/(E/v_phi) are 1.16318, 1.28186 and 1.37837 for the
three charge states. This does not invalidate the accepted mass/depletion
identity or the D Born relation. It prevents promoting their parameter
derivative into an unspecified universal fifth force. The proposed portal
avoids that ambiguity by coupling directly to phi^2-v_phi^2.

## Numerical signal-size examples

These examples reuse the classical Q=39,944 profile and unchanged
dimensionless couplings. They hypothesize other scales; they are not the
original frozen benchmark or dark-matter fits. The halo is supplied and
assumed entirely composed of this one prepared charge. The ideal xenon-131
count assumes perfect efficiency, no background and a 5-40 keV window.

| Hypothetical v_phi | Object mass | Energy R90 | Ideal events/(tonne year), kappa=1e-10 | Formal kappa for one ideal event/(tonne year) |
|---|---:|---:|---:|---:|
| 1 GeV | 3.639 TeV | 10.07 fm | .01114 | 9.48e-10 |
| 10 GeV | 36.39 TeV | 1.007 fm | 7.00e-5 | 1.19e-8 |
| 1 TeV | 3.639 PeV | .01007 fm | 7.18e-11 | 1.18e-5 |
| Original v_phi/M=1e-5 | 8.862e16 GeV | 4.136e-13 fm | 4.97e-42 | 4.49e10; outside validity |

The last column is a normalization diagnostic, **not a detection forecast,
an allowed coupling or a GRUT-selected value**. It shows why the original
high-scale point is not rescued by a weak portal. Its formal required
coupling is vastly outside perturbative/stable-quartic applicability.

The low-scale example has an additional shape consequence. For Q=39,944
and v_phi=1 GeV, the object's squared form factor is:

| Xenon recoil | 5 keV | 10 keV | 30 keV | 40 keV |
|---|---:|---:|---:|---:|
| F_Q(q)^2 | .49170 | .23715 | .007706 | .0004443 |

Over the window its ideal count is .15583 of the point-object count at the
same zero-momentum cross section. The high-energy recoil suppression is
therefore a potential conditional shape test. It cannot be inferred from a
standard point-WIMP limit by comparing cross-section numbers alone. With
v_phi=10 GeV the window ratio is .97635, so that example has little such
shape discrimination.

For three prepared charge states, the mass-normalized zero-momentum
cross section is common at the same v and portal coupling, as derived.
An unspecified charge distribution would average distinct recoil spectra;
it cannot be replaced by this monodisperse illustration without an assumption.

## What must happen before this counts as viable empirical physics

1. **Existing-data confrontation.** Fold the extended-object spectrum through
   an actual detector's efficiency, energy resolution, isotope mixture,
   selection and background likelihood. Assess nuclear and halo uncertainty.
   The public 2024 LZ result provides a concrete dated reference, not a fit
   performed here; the collaboration's 2026 analysis is not reinterpreted.
2. **A visible-interaction-consistent history.** The portal also mixes h and
   s, permits new SM production and can populate/heat the hidden sector.
   Where h-to-ss is open it must be included. Recheck symmetry restoration,
   walls, radiation, long-lived scalar decays and relic populations. The
   gravity-only cold histories cannot simply be carried over.
3. **A prepared population.** Formation, surviving charge distribution and
   halo fraction remain to be obtained from a stated history, not from the
   measured dark-matter density alone. The table supplies all three.
4. **The complete matched EFT.** Add the allowed counterterms, specify the
   portal matching, run induced H-X/H-D interactions and rematch SM masses.
   Check collider, flavor, fifth-force and cosmological constraints in the
   appropriate mediator-mass range. The tree-level calculation is not a
   radiatively complete exclusion analysis.

None of these requirements selects the absolute vacuum density. That stays
a separate unresolved obligation. Adopting a calibrated total vacuum term
would be an honest phenomenological input, not a derived solution; no such
change is made here.

## Decision and next bounded calculation

The empirical route is now specified enough to evaluate instead of just
requesting a generic "visible coupling." **The next discriminating task is
a detector-likelihood recast coupled to a portal-heating/decay history, not
another internal identity.** Start with the explicitly hypothetical low-
scale points above, and reject any signal-size interval that fails either.
Do not label one viable before both checks.

Whether to adopt that new portal and scale as a GRUT candidate is an owner
choice. The scientific work here is a feasibility proposal. Candidate 1,
all prior scientific packets and canonical GRUT-RAI are unchanged. See
[VERIFICATION.md](VERIFICATION.md) and [REPRODUCE.md](REPRODUCE.md).

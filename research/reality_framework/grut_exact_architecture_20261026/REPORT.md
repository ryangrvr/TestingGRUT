# An explicit theory architecture, with an exact clock-and-record instance

26 October 2026, America/Chicago. Codex / AInstein.
Parent 9877518b65989c79adb742e1ff33d7d4b0f85a7d.

**The concrete working proposal is a constrained quantum universe with a
covariant effective field action, internal clocks and physical record
effects.** [THEORY.md](THEORY.md) specifies it, including the supplied state
and matching inputs. The new work makes part of that proposal executable
and adds an explicit quantum-matching obligation. It does not select a
unique fundamental law from closure.

## What is now explicit

The theory datum lists state space, constraints, action, coefficients,
preparation and physical observation effects. The local field content is
SM + metric + phi,X,Y,Sigma + hidden U(1), with charges 1,3,6. A compact T^3
branch is a stated choice. The scalar potential basis is fully enumerated;
allowed quantum counterterms are distinguished from zeros chosen at one
matching scale. Vacuum density is a renormalized input. Quantum probabilities
use the Born rule. Microscopic gravity and the cosmological preparation law
remain unconstructed.

Two calculations accompany that specification.

**First, an exact closed clock/record model.** A 64-dimensional kinematic
space contains an eight-reading clock, one matter bit and two record bits.
The one time-independent constraint has an eight-dimensional physical
kernel. The physical effects commute with the constraint and give the
correct conditional probabilities. At readings 3–5 two records contain the
same one bit of pointer information. A whole-system interference witness
still distinguishes the coherent state from a mixture. Later the records
are erased. The source circuit and blank input are supplied.

The same constraint permits any record probability p. A maximally mixed
physical preparation produces no changing conditional records, and global
purity leaves all p available. An added explicit state filter selects one
example, demonstrating precisely what additional state selection would do.
This is a worked dictionary and a scoped selection counterexample, not a
new quantum interpretation or a solved observer problem in gravity.

**Second, a quantum operator audit.** There are 26 hidden-scalar invariant
monomials through dimension four, or 19 with exact phi parity, before
redefinitions. Exact rational differentiation gives

\[
\left.\mathrm{Tr}[(M_s^2)^2]\right|_{\kappa^2}
 =36\kappa^2|X|^4+72\kappa^2|X|^2|Y|^2.
\]

The mixed X–Y quartic must therefore enter the quantum matching basis. Its
isolated full-theory beta contribution at the supplied kappa is 2.72712e-6.
The earlier gauge construction already named the missing portal/hierarchy
matching; this result quantifies one piece. It neither changes nor excludes
the pinned low-energy target used by the kinetic solver.

The Landau-gauge MS-bar vector-potential derivative at mu_R=m_A gives one
contribution −e^2 m_A^2/(16 pi^2) to the coefficient of |X|^2. At fixed h
its magnitude grows as e^4. This demonstrates why lowering the vector relic
by raising e still requires a stated light-mass matching. The contribution
is not a physical pole correction or an imposed naturalness criterion.

## A concrete covariant replacement for the prescribed Hubble mass

One specified representative is xi_X R|X|^2 with xi_X=−5/6. It reproduces
−10 H^2 during ideal de Sitter inflation and vanishes in ideal radiation.
In an ideal matter era it instead gives −2.5 H^2. This ties that mass term
to the internal metric evolution and fixes its behavior across equations
of state, conditional on adopting this extra coupling.

At the prepared benchmark amplitude F stays positive and the fractional
Planck-mass-squared change is about 2.92e-11. No reheating or phase-holder
trajectory is solved. Replacing the old phenomenological term by this
covariant representative would define a new history to test; it does not
retroactively certify the old finite-release runs.

## What each part explains at this boundary

| Question | Concrete construction or equation | Remaining input |
|---|---|---|
| What exists? | Stated matter fields, metric, gauge group and scalar operator basis | Why these fields, coefficients and topology |
| What changes? | Adopted covariant EFT action; a solved finite propagation constraint | A complete microscopic gravitational dynamics |
| What is time? | Exact conditional observations of the internal finite clock | Physical cosmological clock/dressing and its accuracy |
| What is a record? | Two internal pointer systems with a checked physical POVM | Macroscopic stability and the state that produces a record arrow |
| What is matter organization? | The prior conditional scalar bound-object construction | Formation, full lifetime/population and observed evidence |
| What sets the cosmic history? | A preparation datum and an explicit curvature-mass option | Physical phase holder, reheating, state selection and the kinetic result |
| What is dark energy? | Independent renormalized vacuum stress in the Einstein equation | A selection law for its value or evolution |

The outstanding GRUT claim has a precise location: a restriction on the
action, preparation or physical probability rule with an operational
consequence. The finite state filter illustrates strictness; its choice
is not earned by calling all baths internal. Renaming the existing rules
would not supply that missing step.

## Verification actually run

Eight mathematical groups pass: four for the exact finite model, and four
for operator completeness, rational Hessian, vector-potential derivative
and the curvature-mass representative. The exact spectrum agrees within
1.78e-15; physical conditional probabilities within 1.45e-15; effect
constraint commutators within 3.21e-16. Forty scalar Hessian samples agree
with independent gradient differences within 2.84e-10. The gauge derivative
finite-difference error is 2.90e-11. These are implementation checks,
separate from physical viability or empirical evidence.

No new production release, kinetic solver, inflation scan, detector fit,
macroscopic-record simulation or quantum-gravity quantization is run.
Claude retains the kinetic calculation. See [REPRODUCE.md](REPRODUCE.md),
[HISTORY_RESULTS.json](HISTORY_RESULTS.json), [OPERATOR_RESULTS.json](OPERATOR_RESULTS.json),
[ACTION.md](ACTION.md), [RECORDS.md](RECORDS.md) and [SOURCES.md](SOURCES.md).

# Verification actually performed

One author, Codex/AInstein. No independent reviewer, subagent, canonical
merge, discovery claim or frozen predictive test. Old 51/99-check packets
were read as inputs; their entire checks were not rerun here.

## Experimental implementation controls

The supplied 50 GeV SI spectrum at 1e-45 cm^2 gives 173.560165 accepted
events in the original public implementation. Our 16-point integration
gives 173.560163. Its asymptotic limit is 5.36742782e-47 cm^2, versus
our 5.36742785e-47. The recommended-threshold limits differ by 0.0163%,
within the declared 0.05% implementation tolerance, and far below the
public likelihood approximation's physical/statistical accuracy.

The unmodified official implementation also computes the reference
object's observed limit and the first three future toys. Its future
strength limits are .477326195487111, .226767296348345,
.157014104564716; our vectorized combination agrees below 1e-9
relative. All 1000 toys are used for each of the twelve projected cases.
The largest observed strength-limit change from 1000 to 4000 grid
points is .0301%; 8-to-16-point recoil-bin integration changes a
maximum-bin-normalized value by at most 2.03e-12.

The migration matrix normalization was explicitly checked. An early
development run omitted the future factor of 20 in exposure. That was
corrected before the archived outputs; the public-code comparison
checks the factor using exactly the same reference spectrum.

## Production and history controls

The boosted two-body inverse-energy moment is integrated directly in
rapidity, separately from the Bessel-K0 hot source implementation;
relative residual 2.22e-16. Direct spin-averaged ee->XX* phase space
and the width/spectral-function route agree to 2.22e-16. These are
numerical formula controls, not evidence for the action.

The cold integral changes by 4.00e-7 relative when temperature samples
increase 401->801 and invariant-mass samples 5001->10001. The explicit
2-1000 GeV tail under Gamma_SM<=m is 3.14e-48 GeV per entropy,
negligible compared with the source uncertainty. DOP853 and Radau hot
integrations agree on peak inferred equilibrium temperature to
3.86e-11 relative. Conserved comoving radiation energy closes to
2.22e-16 for the reference. None of these controls establishes rapid
hidden chemical equilibrium or computes an object-formation history.

## Original nuclear implementation

Ran three original ACROPOLIS 1.3.1 decay calculations at effective
120 energy / 20 temperature points per decade: nominal source,
enlarged source and an additional long-lived pure-photon case. The
last is rerun at effective 150 / 30 resolution. The change in the
computed D/H increment is .2815% and in the He3/D increment .2809%.
Those are relative errors on already small increments. All runs stay
inside the package's supplied cosmological temperature range.

During development, changing only `params.NE_pd/NT_pd` after imports
left module-local copied settings unchanged. The wrapper was corrected
to set `cascade.NE_pd` and `nucl.NT_pd` too; the archived finer run
uses the actual 150 / 30 grids and has an explicit settings assertion.
The other three runs already used the package defaults 120 / 20.
No installed ACROPOLIS source file was edited. The runtime entries
are non-scientific timing metadata.

These runs compute incremental photodisintegration on the supplied
standard background; they do not refit initial BBN, CMB or latest
nuclear data. Scalar populations are energy-budget overestimates,
with the branching/lifetime stress choices explicitly disclosed.

## Fresh isolated replay

`verify.py --replay` copied the calculator, compact public arrays and
specified parent profiles to a temporary tree and freshly executed the
detector, history and shape calculations. It reproduced 213 detector,
699 history and 17 shape scientific numerical leaves within 3e-10
relative / 1e-35 absolute tolerance, five CSV tables and the numerical
arrays in both NPZ archives. The archives are compared by array values,
not container timestamps. This includes all 13000 entries of the
1000-toy coupling table. No duplicate full nuclear replay is claimed;
the separate original nuclear runs and effective-grid refinement are
the nuclear verification. A later prose-only scope clarification in
the history calculator/results did not change numerical values.

`VERIFICATION_RESULTS.json` records the original-code crosschecks,
distinct formula routes, tail integral and environment versions;
`REPLAY_RESULTS.json` records the fresh replay scope. Full source/data
hashes, Python compilation, finite JSON/array values, Markdown local
links, figure rendering and Git whitespace are checked before release.
The final scientific figure was rendered and visually inspected.

Installed dependencies were isolated in a temporary target directory:
ACROPOLIS, numba, llvmlite and tqdm. The existing runtime supplied
NumPy/SciPy/Matplotlib. This was not a no-install calculation. See
[REPRODUCE.md](REPRODUCE.md) for exact versions and commands.

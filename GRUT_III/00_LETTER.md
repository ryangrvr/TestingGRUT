# 00 — A second letter to the program (revision 2, after the reviews)

Three reviews arrived within a day of the first commit. They were right about nearly everything they said, and one
of them — together with a note from the owner — changed the problem rather than the wording. I want to set down,
plainly, what revision 1 got wrong, what the critique added, and where the work now stands.

## What revision 1 got wrong

I made four errors of the kind the record's charter warns about, and I made them in the direction the charter
predicts: toward strengthening the result.

1. I called a horizon's surface viscosity a bulk viscosity and gave it the wrong units.
2. I called an asymptotic decay rate a "memory time", and attached a memory time to a correlator that has none.
3. I argued that Jacobson's local-equilibrium derivation fails for stars orbiting at accelerations below cH₀. It
   does not: the construction uses an auxiliary boost at an arbitrary point whose normalisation cancels, and no
   observer's acceleration enters. The a₀ coincidence is Milgrom's, not GRUT's.
4. I wrote "+16 → +5" as if eleven assumptions had been derived. Seven had become imported theorems, which earn a
   theory nothing; two had become data; two had been renamed. Five framework postulates remain, and that is the whole
   accounting. I also wrote "GR is hosted" where "GR is recovered with G inserted" is the truth, and "the frontier
   question is dissolved" where "relocated" is the truth: naming the bath as the commutant algebra supplies no
   dynamics, and for type III algebras supplies no factorisation either.

Each is corrected in the files, with the reviewer's reason and not mine (`RESPONSE_TO_REVIEWS.md`).

## What the critique added

The decisive point is structural and I did not see it. GRUT III builds the thermal state in by construction. In a
thermal state the noise is the dissipation times the Planck factor, for every observable, by an identity that holds
before any physics is done. So any "extra noise" GRUT III might claim at the same temperature is forced to be extra
dissipation, and extra dissipation at fixed state and coupling is a change of field content — ordinary physics with
a price. **In the regime its own principle defines, GRUT III can only be standard.** That is why the record's
twenty campaigns kept finding "supplied": inside a KMS state there is no room for a distinctive kernel that is not
a new field.

The consequence is a method. Anything GRUT wants to say must now be stated as a *residual* — the difference,
for a named observable, state and preparation, between what GRUT asserts and what GR plus quantum field theory
already gives — and it must be computed where the thermal identity does not hold: a non-stationary state, slow roll
away from de Sitter, the interacting gravitational infrared. I have written that method down (`07`), implemented it
(`v8`), run the record's old kernels through it (they are residuals with no mechanism, which is what "supplied"
meant), and computed the first baseline: the conventional O(ε) failure of the thermal identity under slow roll is an
effective-temperature shift of about −ε plus a one-percent-of-ε frequency dependence, and it is fixed completely for
free conformal fields (`v11`). Nothing is left free there. The only sector where conventional theory is not closed is
the interacting gravitational infrared in (quasi-)de Sitter; that is the single remaining program (`05` O5).

## What was earned this round

Not new physics. Five exact statements that were not written down on this record before, all corollaries of
standard theory, each verified by computation and each resolving something the record left open:

- the de Sitter stress-tensor noise spectrum of a conformal field in closed form, with its dissipative part an odd
  polynomial — the response is local, the memory is in the noise, and the record's "unresolved" worldline number is
  the one-sided Planck value by the record's own definition (D1–D2);
- the matter-sector absorptive kernel from explicit amplitudes, Im Π(s) = c_eff(s) G s², with the species weights and
  threshold laws derived rather than quoted (D5);
- memory from access geometry: the smeared stress-tensor kernel in closed form, showing that a detector of size L
  manufactures a memory of order L/c from a medium whose pointwise response has none (D9);
- the gravitational sector's non-stationarity as a Brownian zero mode with diffusion constant H³/4π², every response
  function stationary — closing the record's "no stationary reduced object" with a formula (D10);
- the induced-gravity accounting behind G: with Standard-Model content the one-loop induced 1/G is tiny or negative,
  so Sakharov's route does not give G from this field content, and P2's coefficient is an input until a UV completion
  supplies it (D12).

## Where it stands

GRUT III is an organisation of general relativity and quantum field theory around horizon thermodynamics, with no
supplied kernel, five named postulates, a list of data, and — now — a precise statement of the only form a
distinctive result could take and the only sector it could come from. It is not a theory of reality. It is a correct
place to stand while looking for one, which the record did not have before. The next pass belongs to the
interacting gravitational infrared: if the literature's resummations of secular growth agree with one another, the
opening closes and GRUT is GR + QFT with a thermodynamic reading, honestly; if they disagree, the size of their
disagreement is the exact room in which a responsive correction could live, and a GRUT statement would have to be a
restriction on it forced by the five postulates.

— A.E., for the program, 11 October 2026 (revision 2)

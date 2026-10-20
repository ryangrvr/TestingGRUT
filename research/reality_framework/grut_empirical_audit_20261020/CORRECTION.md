# Qualification of the earlier empirical-exclusion report

The report at commit `479f111` correctly rederived the leading charge-energy
floor and its conditional flat-FRW H0 bound. The same-author audit reproduces
its nominal and conservative numerical bounds to roundoff. The benchmark is
retained as excluded within its specified physical and observational scope.

Its description of the megamaser H0 comparison as independent of cosmology or
the CMB was incomplete. Pesce et al. 2020 Eq. 1 converts geometric disk distances
and redshifts into H0 using flat LambdaCDM and Omega_m=.315 from Planck. The
abstract's independence from the CMB/distance ladder refers to the distance and
velocity measurement basis; it should not be expanded into an exactly
cosmology-independent inferred H0. The old algebraic H0 floor is independent
of hidden pressure; this does not make the observed H0 estimator independent
of its assumed expansion history.

The public freeze `eb5a1084cb8e9eaf4bef8c6cfa1f9cbd8b93d153` specifies a
replacement **comparison observable**, without changing the action, mass or
yield: compare the conservative D_A(z) envelope directly with all six Table 1
disk distances. This avoids the Planck Omega_m prior and a hidden-dust EOS.
It still assumes flat weak-binding FRW, standard BBN and the stated disk and
peculiar-velocity interpretation. It is not free of every cosmological input:
the velocities use the CMB rest frame, and the baryon-number conversion adopts
the photon-temperature/standard-BBN history. No CMB anisotropy density prior
is used in the new prediction.

The distance comparison uses the **same maser observations**, not a second
independent test added to the old H0 statistic. The table and old rejection were
already known before the source freeze. No blinded, prospective, independently
replicated or positive-identification claim follows. Conservative percentile
corners are sensitivity intervals, not a calibrated joint confidence set.

One further implementation clarification: the mathematical `prediction()`
function reads protocol and calibration values only. The command first runs an
integrity guard that hashes all frozen files, including the test-table bytes.
It is not an information-blinding mechanism, and the table was already exposed.
The estimator itself consumes no test-distance, velocity or fitted H0 value.

The earlier frozen packet is left byte-identical so its original statements,
failed controls and numerical records remain reviewable. Read this correction
and the new `REPORT.md` with that packet. No old data are erased or refitted.

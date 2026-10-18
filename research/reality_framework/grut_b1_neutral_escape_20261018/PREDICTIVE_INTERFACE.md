# B1 observable and calibration interface

This is a proposed interface for an engineered quantum medium, not an
experimental result. No laboratory data were used, no device is claimed to
exist, and no prediction distinguishes GRUT from ordinary quantum mechanics.

## Dictionary and units

| Model quantity | Possible measured quantity | Supplied mapping |
|---|---|---|
| Fuel F | Population of an initially excited neutral register | Register implementation |
| Plus/minus occupancy | Two distinguishable hard-core excitations | Charge-sector encoding |
| E on a link | Three-level electric register | Synthetic constraint implementation |
| Co-located pair | Coincident opposite-charge occupations | Joint occupation readout |
| Recorded distance x | Number of links along the recorded word | Relational distance, not meters |
| Color | Internal two-level degree on fuel or plus | Coherent encoding and decoder |
| theta=Jt | Calibrated evolution time | Hamiltonian energy scale |

If E_* is the energy multiplying the dimensionless Hamiltonian, physical
time is t*hbar/E_*. Thus preparation time is 7680*hbar/E_* and final readout
15360*hbar/E_*. A physical length per recorded link requires an additional
device mapping; it is not a fundamental length predicted by GRUT.

## Calibrate, freeze, then predict

Calibrate constituent hoppings, mass, electric coefficient, density coupling,
fuel frequency/coupling and endpoint counterterm using separate local
spectroscopy. Fix which record, source, time units and arrival set are used.
Reject an implementation outside the declared parameter class rather than
fitting an arrival trace to move it back inside. Parameters not covered by
the tolerance class must be known to the precision required by an additional
error account; the ideal theorem does not supply a generic device-noise model.

Do not use arrival probability or retention at theta in [2,4] as calibration
data. These are the holdout outcomes:

* p_matter>1/2 throughout [2,4]. A finite time sampling alone cannot verify
  the entire interval without a separately bounded interpolation error.
* At theta=4, preserved-prefix co-pair arrival at x>=4 succeeds with
  probability >3/5 and preserves arbitrary input color on success. All other
  trials are recorded erasures.
* The nominal reference fuel amplitude is J_1(2theta)/theta. Reference
  arrival profile is (n+1)^2 J_(n+1)(2theta)^2/theta^2, n=x+1.
  For a projection in the fuel/co-pair subspace, its physical probability
  lies between max(0,sqrt(p_ref)-e)^2 and min(1,sqrt(p_ref)+e)^2,
  with e=e0+rho*t+e_end+2theta/1000 plus the full-geometry D_B if applicable.

One can compare return probabilities, co-pair arrival distribution, erasure
rate and color tomography. Statistical confidence and readout error would
have to be added before an experimental decision. The theorem's approximation
error is not a substitute for those experimental uncertainties.

## Conventional comparison

The reference is a standard continuous-time tight-binding walk with one
endpoint initially occupied. It predicts the same Bessel profiles. The
intertwiner accounts for virtual charged excitations in B1, not a correction
to quantum mechanics. A successful fixed-path experiment would support that
engineered Hamiltonian and approximation, not GRUT's geometry-formation
ontology. A failure inside a genuinely calibrated class would reject the
claimed implementation/model correspondence or its derived bound.

The coherent one-composite result is a candidate mathematical/quantum-medium
prediction protocol. A distinctive prediction about nature remains none.

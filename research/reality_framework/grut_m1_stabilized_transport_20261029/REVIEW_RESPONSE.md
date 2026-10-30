# Response to the supplied Argon/Gemini proposal

The proposal usefully identified three calculable extensions. It was not independent numerical replication. This response records what the new calculations establish and which suggested claims fail.

| Proposal | Checked outcome |
|---|---|
| Replace the sextic curvature holder with `R Sigma* Y²/M_Y` | Implemented. Full mixed-quartic bound is positive at leading scalar degree; explicit finite-R lower bound proved. |
| Cubic scalar degree establishes UV stability/renormalizability | Does not follow. The curvature operator has dimension five; scalar boundedness is not nonnegative energy or global gravitational/quantum consistency. |
| Prescribed phase becomes a state-selection law | The selected complex coefficient encodes the phase. A physical initial-state mechanism remains missing. |
| `|R||X|¹²` term is an exact square completion | Not a square with the Sigma quartic. It can still stabilize the fixed-R scalar potential as explained in DERIVATIONS. `|R|` is nonanalytic at zero curvature. |
| `R²|Sigma|²|X|⁶/Lambda¹⁰` stabilizes the sextic potential | Wrong operator dimension for dimensionless c8; correcting to Lambda⁸ still leaves a fixed-Sigma unbounded X direction. |
| Reuse the 12×12 engine for actual CMB isocurvature | Its scalar/longitudinal transfer was extended, but the original engine excludes metric, inflaton and radiation perturbations. Computed scalar charge power is conditional, not actual CDI. |
| Estimate decoupling near 0.62 MeV | Homogeneous viscosity-proxy crossing is 0.62156 MeV in the implemented model. This is not a universal perturbation decoupling surface. |
| Eight to sixteen collisions ensure thermalization | Approximately 7.91 collisions are computed, and the final shape is measured to be close to Maxwellian. No universal collision-count theorem was assumed. |
| Smooth the final PSD to make CLASS converge | Some moderate-k Runge–Kutta controls succeed. A high-k warm NDF run still fails, and the cold comparison times out. |
| Produce a parsec-scale transfer function | A sampled conditional BGK test-cosmology curve is now calculated and refined. Full interacting M1 power remains open because the formation perturbations and exact collision hierarchy are missing. |

Correction to our previous report: the quoted 0.841% hypothetical mass correction need not be an arithmetic mistake. It corresponds to cooling the existing gas heat completely, whereas 0.85755% retains that heat. Both additionally assume the unproved conversion of the remaining reservoir. This packet keeps both histories conditional and does not rewrite the old packet. Its abundance table uses the fully cooled 0.84099365% hypothesis where explicitly labeled.

Verification history is retained honestly: the first replay criterion `1e−8` failed on a `1.96e−8` difference after extraction of the algebraically equivalent collision routine. A separately tighter replay differed by `9.42e−9`, while combined grid/hierarchy/radiation averaging shifts are about `3e−5`. The replay criterion was revised to `1e−7` with that reason recorded in VERIFICATION.json. No global `1e−8` power accuracy is claimed.

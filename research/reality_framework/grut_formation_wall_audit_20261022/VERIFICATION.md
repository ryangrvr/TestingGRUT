# Verification and limits

One author, no independent replication or subagents. No new packages
were installed for this packet. Old detector, nuclear and wall-network
runs are inputs/context; they are not represented as rerun here.

The following numerical controls were actually executed:

* Direct kink energy integration agrees with the analytic neutral
  tension to machine precision. Direct Yukawa trace integration agrees
  with 2 mu alpha/m_s to 2.22e-16 relative.
* An independent organization of the full radial/phase/phi linear
  equations as a 6-by-6 first-order matrix checks the growth determinant
  at twelve denser control configurations. These avoid trying to
  subtract O(m_X) eigenvalues to extract a 1e-13 GeV physical rate.
  Stable outside-band controls have no positive real growth above
  the declared numerical tolerance. The largest relative growing-mode
  difference is 2.74e-11.
* A separate Radau integration evolves the full linear equations for
  five growth e-folds, with both massive fields and phase retained.
  Its growth agrees with exp(5) to 4.93e-12 relative. This is a linear
  control, not a nonlinear evolution of formation.
* At the actual reference density, a 10001-point scan around the
  maximum reproduces the optimized peak. The small-k sound speed
  approaches the independently derived charge-compressibility result
  with differences 5e-4, 5e-6 and 5e-8 as k^2/(a f_0^2) is reduced
  from 1e-3 to 1e-5 and 1e-7. This is the expected finite-k correction,
  not an estimate of the model's physical error.
* Charge/energy accounting is checked across all twelve previously
  supplied object/scale cases. Reusing them does not derive their
  charge distribution, calibrated abundance or allowed status.

`verify.py --replay` copies only the calculator and its two frozen
parent JSON inputs into a temporary isolated tree. Fresh execution
reproduces all 522 scientific numeric leaves and both CSV arrays at
3e-10 relative / 1e-35 absolute tolerance. This is same-author
reproducibility of the declared calculations, not 510 independent
physical tests.

The manifest pins delivered files and both consumed parent inputs.
Before publication, Python compilation, finite JSON/CSV values, local
Markdown links, hashes, Git whitespace and the rendered scientific
figure are checked. Published Git blobs/tree and report readback are
checked separately; no canonical merge is made.

Not calculated: primordial charge production, gas-to-condensate
preparation, nonlinear fragmentation or a final mass function, a
formation radiation spectrum/branching ratio, a new nuclear or CMB
likelihood, a quantum lifetime, new wall-network evolution, or a
multi-target detector forecast. The total vacuum matching remains a
calibrated assumption. No result here supplies today's vacuum density.

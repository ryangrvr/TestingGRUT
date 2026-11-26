# F1 frozen analytical protocol

26 November 2026. Parent is draft PR #23 at
14cf11c7fed1413d63f9139f25be9e1569a97f8c. The following proposition,
controls, observational interface and exact calculation budget are fixed
by this file and WITNESSES.json before the executable checks.

## Proposition and domain

Take a connected, finite, undirected spring network with positive
edge weights, labelled canonical site displacements and momenta,
unit masses, and a known common pinning frequency squared of 1.
The stiffness is $V=I+L$ for the physical spring Laplacian $L$.
The Hamiltonian is $H=(p^{\mathsf T}p+x^{\mathsf T}Vx)/2$ in
the declared units. A nonempty set of calibrated sites can be
driven and measured. Let $P$ select those sites. The allowed data
are the entire matrix of retarded displacement susceptibilities
$\chi_P(t)$ for positive physical time. The source Hamiltonian
is $-f(t)^{\mathsf T}P x$, so both force and time are calibrated.
The primary proposition assumes that the Krylov span of
$(V,P^{\mathsf T})$ equals the whole site space. This is the
condition that no normal mode is invisible to the data.

**Target.** Determine whether the unknown weighted graph is
identifiable up to permutations of the inaccessible sites. In
this class exact response is expected to identify the smallest
dynamical realization up to an orthogonal change of coordinates.
Its geometry should be identifiable with at most two unmeasured
sites. Test whether three unmeasured sites already allow exact
nonisomorphism, even when all modes are visible, and whether extra
measured sites can remove that ambiguity in general. The
two-dimensional and three-dimensional terms here count unmeasured
site coordinates. No ambient spatial dimension is assigned.

## Matched witnesses and controls

* The cycle and pendant triangle in WITNESSES.json have integer
  springs, one measured site, the same four modes, equal pinning
  and mass, and equal edge count and total spring budget. Verify
  the listed orthogonal transformation, stiffness spectra and
  equality of the entire port matrix. Confirm nonisomorphism of
  the supports, including different triangle counts. Record the
  Krylov rank and the primary absolute determinant.
* Append the registered chains of observed sites and the
  independently attached extra hidden leaf, using the same
  couplings in both networks. Check exact equality of the whole
  expanded response matrix and full rank in each case.
* The connected control with two equal hidden leaves must still
  match response and must show a dark mode, with rank five of
  six. Its failed full-rank certificate stays in the record.
* Check the hidden relabeling control. A permutation fixing the
  old port gives identical data and the same graph. The
  registered calibrated added probe at site 3 must distinguish
  the pair, with third-derivative difference 30 and a certified
  gap greater than $29/50000$ at $t=1/20$. This control raises
  access while keeping the microscopic dynamics fixed.
* For a continuum scope check use the conformally coupled free
  massless scalar in four spacetime dimensions. Compare
  Minkowski with a smooth conformal deformation equal to one
  throughout the accessible open set. Derive equality of the
  restricted retarded kernels with unchanged calibration. This
  is a standard conformal covariance check on prescribed metrics.

Full thermal covariance and any interventions restricted to the same
measured oscillator algebra may be considered under unitary-related
preparations. No extra access to a hidden local state or to an
independently calibrated hidden measurement is presumed.

## Method and decision

Freeze the sources and class before code. Supply a proof of the
positive classification, a connected exact counterexample to
unrestricted identification, and an independent internal critique.
Check the finite witness with integer or rational arithmetic, with
no floating-point acceptance tolerance. One implementation verifies
conjugacy, response moments through order $2n-1$, rank and the
registered separation bound. An independent implementation checks
the projected resolvents and the rank without importing the first.
Each solver gets 30 seconds and at most eight sites. There is no
parameter search. Preserve every negative control.

Verdict C requires a rigorous counterexample within the declared
class and a passing exact replay. A separate result earns A only
if it survives comparators in realization theory, inverse spectral
geometry, causal structure and covariance. A positive conditional
reconstruction with established methods can earn B. D applies if
the stated hypothesis remains underdefined or the claimed
obstruction fails. The literature gate determines the kind of
production permitted. An exact verification of the obstruction
is within budget. A simulator-driven emergence campaign requires
a distinctive mechanism surviving the gate.

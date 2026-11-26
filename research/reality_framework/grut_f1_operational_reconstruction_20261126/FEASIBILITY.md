# Analytical feasibility and stopping decision

The operational hypothesis is that calibrated physical responses can
determine geometry. The test gives it unusually strong data. The
response functions are exact at every time, the forces and site
operators are identified on the measured set, the standard spring
law and common pinning are known, the full set of modes couples to
those measurements, and one compares matched connected networks.

## Prior accounts of this inverse problem

The input-output realization theorem determines a minimal linear
model up to similarity. Symmetric second-order dynamics with a
reciprocal source and readout admit an orthogonal version, obtained
from the Gram matrix of the response-generated Krylov space.
The inferred basis is a mode description. Coordinate changes that
fix all ports can persist on the inaccessible sites.
[Ho and Kalman](https://ntrs.nasa.gov/citations/19670049337)
· [Gough and Zhang](https://doi.org/10.1016/j.automatica.2015.06.023)

Inverse graph spectral reconstruction requires structural assumptions.
For example, Blåsten, Isozaki, Lassas and Lu identify conditions on
boundary adjacency and extreme vertices supporting uniqueness.
Their separate heat-equation result connects no invisible
eigenfunctions to equivalence of two data formulations. The
geometric uniqueness theorem uses more than that condition.
Orthogonal switching is also a standard source of isospectral graphs.
[BIL, inverse graph Laplacian](https://doi.org/10.4171/jst/455)
· [BIL, heat equations](https://doi.org/10.1137/21M1439936)
· [Godsil and McKay](https://users.cecs.anu.edu.au/~bdm/papers/GodsilMcKayCospectral.pdf)

GRUT G-2 already shows a single-site realization ambiguity using a
Lanczos chain. GS-1 tests hidden rotations, reports dynamic separation
of a static Y–Δ degeneracy, and corrects an apparent passivity
uniqueness that depended on the original network. F1 narrows that
existing obstruction to equal masses and **exactly matched common
pinning**, matches the number of positive springs and total
strength, and checks what happens as the observed fraction grows.
That sharpening concerns the inverse problem for a particular
class of ordinary Hamiltonians.
[G-2](../../../G2_SPECTRAL_VERDICT_01.md)
· [GS-1](../../../GS1_GEOMETRY_SELECTION_VERDICT_01.md)

## Why the bounded exact check can decide the question

The calibrated matrix susceptibility generates the exact moments
$P V^j P^{\mathsf T}$. Their joint Gram matrix reconstructs the
reciprocal minimal model. Because $L=V-I$ is a connected Laplacian,
any intertwiner that fixes the measured coordinates must also fix
the constant vector. On a hidden block with one or two entries, an
orthogonal matrix fixing the constant vector is a permutation.
A hidden block of size three supports further transformations.
The registered rational reflection fixes the constant vector and
maps the cycle stiffness to the pendant triangle stiffness while
retaining physical spring signs. The no-dark-mode test and evidence
for multiport equality can both be checked with exact small matrices.

In the continuum, a massless conformally coupled scalar's causal
propagator has the established conformal transformation law. If the
metric deformation vanishes on the apparatus region, neither the
calibration nor the restricted propagator resolves it. The geometry
can differ in curvature elsewhere.
[Pinamonti, conformal covariance](https://doi.org/10.1007/s00220-009-0780-x)

## Gate result

The proposed reconstruction principle has an established
realization-theory and conformal-field explanation. F1 therefore
proceeds with the obstruction proof and a small exact certificate.
The gate closes the novelty-directed construction route at the
declared access interface. No candidate change to field equations,
couplings or initial conditions is warranted by this hypothesis.
The result to verify is the sharp boundary within the stated spring
class, together with its limits as a statement about spacetime.

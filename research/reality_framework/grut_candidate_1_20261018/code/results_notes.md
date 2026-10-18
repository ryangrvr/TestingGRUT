Numerical completion notes
==========================

Reproduce the verified scalar boundary-value calculation and all main observables:

    python GRUT_Candidate_1/code/qball.py

The main script uses the saved numerical preparation. If absent, it executes the included `results_preparation_search.py`, whose first predeclared fixed-frequency bag guess converges. Optional supplementary unconstrained-Hessian diagnostic:

    python GRUT_Candidate_1/code/results_additional_checks.py

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. No extra packages installed. All numerical jobs completed; none remain running.

Physical prescription was held fixed throughout: a=1, b=0.1, g=0.1, hence lambda=0.01 and lambda_X=0.001. Gravity is evaluated perturbatively with v/Mpl_reduced=1e-5, Gv^2=1e-10/(8pi). The successful preparation was selected by numerical initial-state search, then frozen before four verification solves: F(0)=1.979050307129842. Six unsuccessful earlier high-central-amplitude starts are fully archived in `results_initial_*` and `results_gaussian_*`. The fixed-frequency preparation search is recorded separately in `results_preparation_search.json`. No physical couplings were searched or tuned.

Final dimensionless scalar result: Omega=0.8000000019763, S(0)=0.0490799788706, Qbar=399.4384624223, Ebar=363.9229985517, Ebar/Qbar=0.9110865196725, R90bar=5.10343221347. The last two domain/tolerance refinements changed Ebar by 5.87e-8 and Qbar by 7.48e-8. The last solve used L=50 and residual tolerance 1e-9.

Stability scope: real amplitude sectors ell=0,1,2,3; imaginary phase sectors ell=0,1,2; three predeclared finite-difference grids. Fixed-Q rank-one correction is included only for real ell=0. At L40,N3200, lowest eigenvalues: amplitude ell0=0.3669665408 (finite-box continuum level; continuum threshold 0.36), ell1=-2.66211e-6 (translation artifact converging quadratically to zero), ell2=0.2854898742, ell3=0.3905134380. Phase ell0=-1.01340e-6 is the U(1) zero-mode discretization artifact and also converges quadratically; ell1=0.2717298479, ell2=0.3807105444. The local amplitude ell0 operator without the fixed-Q rank-one term has a negative eigenvalue -0.18885499, confirming that the constraint is essential. The exact continuum phase operator has zero solution rho*F>0 and admits the usual nonnegative ground-state factorization. Higher angular sectors add a positive centrifugal term. These calculations provide local classical numerical evidence, not a proof of global nonlinear stability, fission resistance, quantum stability, or stability of unspecified extra interactions.

Stress and gravity: dimensionless e=grad+time+U, pr=grad+time-U, pt=-grad+time-U, grad=(Fprime^2+Sprime^2)/2, time=Omega^2*F^2/2. The dimensionless active source is e+pr+2pt=2*Omega^2*F^2-2*U. Potentials solve Laplacian Phi=4pi*Gv^2*active and Laplacian Psi=4pi*Gv^2*e with zero boundary value at spatial infinity. Metric convention is ds^2=-(1+2Phi)dt^2+(1-2Psi)dx^2. Sources use the flat-space scalar profile; this is the leading weak-gravity metric, not a fully backreacted Einstein-scalar BVP solution.

The stress-conservation residual is 4.45e-10 (maximum absolute), virial/E=1.82e-12, and integrated active mass equals energy to 1.82e-12 relative. Phi(0)=-6.68812557e-10, Psi(0)=-5.23911550e-10, center-to-infinity gravitational redshift=6.68812557e-10 to first order. 2GE/R90=5.67462657e-10. At impact parameter rho_b=20, light bending=2.89600720545e-10 radians, agreeing with 4GE/b to 2.63e-10 relative. Exterior potentials coincide up to numerical error. F decays with kappa=sqrt(1-Omega^2)=0.6; the homogeneous S mass is sqrt(2), but the asymptotic S deviation may be driven by F^2 with exponent 1.2.

Main raw records: `results.json`, `results_profile.csv`, `results_profile.npz`, `results_execution.log`; supplemental diagnostics: `results_additional_checks.json`.

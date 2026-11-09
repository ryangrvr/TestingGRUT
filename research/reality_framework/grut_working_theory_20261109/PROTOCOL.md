# W1: constructive working-theory pass

The task is to consolidate an explicit finite operational theory and construct a usable continuum response approximation. This is not an attempt to force the existing lattice/continuum discrepancies to pass. All historical outcomes remain unchanged.

The independent deliverables are:

1. An explicit calibrated-parameter record-action model. Its general probability law is checked against500 full four-qubit density-matrix evaluations, finite ideal pulse generators, controls and uncertainty examples. These theoretical regression results were inspected before the suite freeze; they are not held-out physical evidence.
2. A new1152-dimensional autonomous dilation, using a nine-position clock and an explicit environment ancilla for the declared parity channel. Parent geometry, four parameter cases, four selector programs, four clock times, probability floor and1e-10 tolerance are fixed in parent_protocol.json. No parent evolution was evaluated before the freeze.
3. A source-mode-resolved continuum predictor. Its primary count1024, diagnostic counts16/64/256, new response times, fixed references16384/32768 and separate numerical gates are fixed in spectral_protocol.json. No W1 physical-source response was evaluated before the freeze. SPECTRAL_DERIVATION.txt proves an exact-arithmetic all-time truncation envelope using rational inequalities and unchanged source integrals.

The exact1024-mode response is a specified mathematical predictor. Its uniform-in-time truncation theorem does not assert that ordinary floating-point frequencies remain accurate at arbitrarily long times. A frequency error obeys |exp(-i omega t)-exp(-i omega_hat t)|<=|t| |omega-omega_hat|; physical parameter uncertainty has the same need for a declared prediction horizon. The executable is checked at the frozen finite times, without an interval-arithmetic roundoff certificate.

Run.py requires clean committed source, hashes every study file and imported implementation, and invokes one serial suite under1GiB Windows Job Object and120-second per-calculation limits,360seconds total, with one BLAS thread. It refuses overwrite, preserves failures and performs no new large-lattice solve. The resource limiter is the previously checked N22 implementation.

    python run.py --output NEW_EXTERNAL_DIRECTORY

The physical source dictionary, source offsets, boundary constants and massless continuum frequencies are unchanged. The spectral repair retains more actual dynamical modes; it does not improve the R1 lattice-to-continuum comparison or validate nonlinear detector backaction. The mathematical parent supplies its clock couplings and preparation as inputs. A physical four-qubit realization still requires the separate M1 calibration and held-out measurement contract. No fundamental selection law, consciousness mechanism or cosmological completion is added.

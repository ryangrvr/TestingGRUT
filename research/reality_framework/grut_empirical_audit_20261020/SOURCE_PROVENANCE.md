# Data and source provenance

Sources were inspected on 20 October 2026. No raw-observation reanalysis,
unavailable Claude output, external review or proprietary detector record was
available. All baseline artifacts are pinned to commit
`479f111117ea7162bd7888ff8f3792b2f6c45fdb`.

| Source | Exact use | Independence and limitations |
|---|---|---|
| [Candidate 1 action](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md), sections 2–7 | Canonical fields, current, potential, GR, matching, parameter status | Project postulates; not new empirical evidence. Read in full. |
| [Baseline packet](../grut_3d_empirical_20261020/REPORT.md) | Report, mathematics, all numerical protocols and empirical amendment, calibration/comparison files, source/raw manifests, development notes and checkpoint | Old 29-run results not rerun. 65 manifest-indexed packet files, three preserved parents, 20 archives and 190 raw members verified. Integrity is not independent replication. |
| [Mossa/LUNA 2020 author manuscript](https://www.pure.ed.ac.uk/ws/portalfiles/portal/166986970/mossa_nature_R2.pdf), DOI 10.1038/s41586-020-2878-4 | Table 1 D+3nu-with-LUNA row and Bayesian-likelihood Methods: omega_b and uncertainty, Neff=3.045 | Primary calibration. The selected row uses deuterium/nuclear data, not the separate D+CMB likelihood. Standard-BBN inference retained. Table visually checked on PDF page 21. |
| [Pesce et al. 2020 v2](https://arxiv.org/pdf/2001.09213), DOI 10.3847/2041-8213/ab75f0 | Table 1 six distance/velocity summaries, Eq. 1 cosmography and section 3.1 peculiar-velocity prescription | Primary test table already seen before freeze. PDF page 4 visually checked. Combined H0 is a contextual result using the same measurements, not additional data. No raw disk posterior samples used. |
| [Fixsen 2009](https://arxiv.org/abs/0911.1955), DOI 10.1088/0004-637X/707/2/916 | Old CMB monopole temperature and uncertainty for SI/yield bookkeeping | Cancels from primary ratio once omega_b is supplied; no new temperature fit. Abstract/value verified. |
| [NIST 2022 CODATA](https://physics.nist.gov/cuu/Constants/Table/allascii.txt) | SI constants in frozen old calibration; conversion cross-check | Declared calibration convention retained; no dark-sector calibration or prediction follows. |
| [Planck 2018 VI](https://arxiv.org/html/1807.06209v4) | Context for old omega_c and the Planck Omega_m used by Pesce Eq. 1 | Not an input to the new distance prediction and not included in a new likelihood. The old cold-matter comparison remains conditional. |

An attempted experimental arXiv HTML route for Pesce returned an access/fetch
error; the primary arXiv PDF was available and used instead. Broad search
results suggesting later maser analyses were not used to substitute data after
the protocol decision. No claim that this historical table is the newest H0
result is made. A fresh independent test would require a reserved dataset and
independently fixed inputs, which are currently absent.

The repository artifacts contain input values, paraphrased assumptions,
equations derived here and source links. They do not redistribute full papers.

A targeted repository search also checked the original calibration/open-claim
documents, the later theory synthesis and physical-access audit for scale,
reheating and yield provenance. It found illustrative scale conversion and
explicitly missing charge/production history, not independent calibration of
the new sector. No parameter rule from an incompatible historical architecture
or inaccessible branch was imported into Candidate 1.

`FREEZE.json` hashes the old calibration and parameter files, the action, the
old raw/source manifests and the original reference-profile results. Every
frozen source file is checked against its public commit by `audit.py`.

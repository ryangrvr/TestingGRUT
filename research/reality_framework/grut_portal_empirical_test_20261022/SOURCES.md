# Primary sources, credit and exact use

This packet is one-author exploratory research, not a literature priority
claim. Citations identify adopted methods and inputs. Computed model
outputs and their limitations are ours; source authors do not endorse
this construction.

1. **XENON Collaboration, approximate nuclear-recoil likelihood.**
   [Paper, arXiv:2210.07231v1](https://arxiv.org/html/2210.07231v1) and
   [official code/data](https://github.com/XENON1T/xenon1t_approximate_nuclear_recoil_likelihood).
   Used release v1.0.0, commit
   `8254684105193d4a84cb8b5752f5e450d5ab70b3`. Read the detector
   response, statistical construction, error discussion and projection
   sections; inspected and executed its original Python implementation.
   The public response is a probability matrix. Its authors describe
   typical approximately 10% heavy-WIMP limit agreement and deviations
   up to approximately 30% for sharper low-energy spectra. The 1000
   XENONnT toy likelihoods belong to the historical 20 tonne-year
   projection, not a contemporary observation. Small public arrays are
   extracted verbatim into `public_inputs.npz`; the 33 MB compressed toy
   ensemble stays external, pinned by hash and commit. No experimental
   code is vendored in this packet.

2. **LZ Collaboration, 4.2 tonne-year search.**
   [arXiv:2410.17036v3](https://arxiv.org/html/2410.17036v3),
   [HEPData record](https://www.hepdata.net/record/155182),
   [SI limits table v2](https://doi.org/10.17182/hepdata.155182.v2/t1),
   [efficiency table v2](https://doi.org/10.17182/hepdata.155182.v2/t4).
   Used the published `limit` column, which is the constrained observed
   point-WIMP limit, and the energy/efficiency table. The retrieved
   JSON endpoints are recorded in `fetch_inputs.py`; their exact bytes
   are hashed. The release was downloaded and parsed, not digitized
   from a plot. Neither LZ calculation in this packet runs the full LZ
   likelihood; both are explicitly count/spectral estimates. Read the
   search description and detector response/limit context. These are
   named dated results, not a claim to have surveyed all data available
   in October 2026. HEPData lists the tabular release under CC0.

3. **Depta, Hufnagel and Schmidt-Hoberg, ACROPOLIS.**
   [arXiv:2011.06518v1](https://arxiv.org/pdf/2011.06518),
   [official repository](https://github.com/hep-mh/acropolis).
   Used installed ACROPOLIS 1.3.1 (GPL-3.0); read the framework/method
   paper and relevant original model, input, cascade, nuclear and
   parameter implementations. Its unmodified cascade/nuclear equations
   and supplied standard-cosmology/BBN data do the disintegration
   calculation. Source hashes are archived. The wrapper changes only
   normal runtime parameters/resolution and turns off verbose output;
   no installed source file is edited or copied into this packet.
   The wrapper sets both parameter-module values and the module-local
   settings actually consumed by the solver. We do not claim an
   independent nuclear code or an updated BBN data analysis. The
   package's supplied observations/background correspond to its
   documented version, not a newly assembled 2026 dataset.

4. **Fradette and Pospelov, Higgs-portal scalar lifetimes and BBN.**
   [arXiv:1706.01920](https://arxiv.org/abs/1706.01920).
   Read the scalar-decay/hadronic-response discussion and abundance
   qualifications. Used established lepton-width and low-energy scalar
   current methods, not its published bounds as universal lifetime
   cutoffs. Our masses avoid physical muon/hadronic s decays; uncertain
   off-shell photon/chiral production is explicitly stressed. C_gamma
   choices 10/100 and hadronic factor 1000 are ours, not fitted values
   or quoted uncertainty intervals from this paper.

5. **Hoferichter, Klos, Menendez and Schwenk, Higgs-nucleon coupling.**
   [arXiv:1708.02245v2](https://arxiv.org/abs/1708.02245v2),
   PRL 119, 181803 (2017). The abstract explicitly gives f_N=.308(18),
   adopted here as a reference. That uncertainty alone changes the
   cross section by about 12% at first order and the inferred coupling
   by about 6%. The detector tables fix f_N and halo parameters; they
   are not nuisance-marginalized confidence regions. No fresh full
   rereading of this paper is claimed.

6. **NIST xenon isotopic composition and atomic masses.**
   [Official table](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Xe).
   Read and transcribed its nine central isotope fractions and masses;
   both are explicit in `model.py`. Nuclear masses subtract the electron
   contribution at the precision used here. The Helm and halo choices
   are declared conventional modeling inputs, not experimental fits.

7. **Barbieri et al., low-reheating constraints.**
   [arXiv:2501.01369v2](https://arxiv.org/abs/2501.01369v2).
   Read its low-reheating definitions and bound context; its abstract
   reports T_RH>5.96 MeV at 95% under its analysis assumptions. This
   motivates investigating a supplied 20 MeV history, but does not
   construct ours. A maximum temperature after preparation and a final
   reheating temperature are distinct conditions. We do not import its
   likelihood or claim that this one inequality validates our cosmology.

8. **Higgs-portal soliton precedent.**
   Ponton, Bai and Jain,
   [Electroweak Symmetric Dark Matter Balls, arXiv:1906.10739v1](https://arxiv.org/abs/1906.10739).
   The parent packet inspected its scattering section and credited it.
   It is a related, different action/scale; its solutions, formation and
   phenomenological survival are not transferred to the present model.
   The Higgs low-energy/background-derivative method, including SVZ
   1978 and the 1980 review, is credited with precise access limitations
   in the [parent sources](../grut_higgs_portal_feasibility_20261022/SOURCES.md).
   No newly read original full text or protected Ward identity is asserted.

The compact data archive preserves selected public numerical inputs,
with attribution, rather than copying article text. `SOURCE_INPUTS.json`
is the byte/array provenance record; `MANIFEST.json` also pins our
parent profiles and each delivered file.

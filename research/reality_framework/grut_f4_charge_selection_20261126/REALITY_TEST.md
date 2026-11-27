# Reality test for electric charge selection

**Observable.** Use the long-wavelength electric monopole of the free $udd$ hadron divided by the magnitude of the charged lepton's electric monopole. Calibrate the electric charge unit separately from the neutron deflection measurement. The test concerns an identifiable long-lived neutron analogue and ordinary electromagnetic coupling. In the chosen comparison class the quark masses and strong coupling can be matched in the QCD-like regime. Interactions so large that the comparable species disappears are outside this empirical comparison.

\[
Q_n(k)=\frac{q_n}{|q_e|}=\frac{3k}{|1+3k|}.
\]
A positive $k$ has $|Q_n(k)|=1-1/(1+3k)\ge3/4$. A negative $k$ has $|Q_n(k)|>1$. The difference from the canonical electromagnetic charge is therefore at least $0.75$ electron charge magnitudes in every noncanonical member. Charge additivity and the Ward identity fix the monopole within the model. Hadron structure can modify form factors beyond the monopole, while the integrated charge is fixed by the representation. A neutral hydrogen atom supplies no bound on $k$ in this family because $q_p+q_e=0$. [Derivations](DERIVATIONS.md)

| Registered $k$ | $Q_n(k)$ | Unbroken-phase bound on $d_*(k)$ | Pass neutron benchmark |
| ---: | ---: | ---: | :---: |
| -2 | $-6/5$ | $23\le d_*\le39$ | No |
| -1 | $-3/2$ | $11\le d_*\le21$ | No |
| 0 | $0$ | $d_*=3$ | Yes |
| 1 | $3/4$ | $15\le d_*\le21$ | No |
| 2 | $6/7$ | $27\le d_*\le39$ | No |

**Empirical input and uncertainty.** Baumann et al. searched for electric-field-correlated deflection of a cold free-neutron beam passing through a 9 m field region. Their Eq. (10) reports $(-0.4\pm1.1)\times10^{-21}$ in electron charge units at 68% confidence. Eq. (9) uses different run weights and gives a center shifted by $0.2\times10^{-21}$ with the same reported uncertainty. The primary comparison uses Eq. (10) as frozen and only the absolute center enters. The experiment tested possible field-correlated optical shifts and reports that uncertainties in wavelength, lengths and applied electric field were small compared with the statistical errors. This packet uses the published combined reported error at face value. It does no raw-data refit and assigns no p-value to a branch many orders of magnitude outside the small-deflection calibration. [Primary full text §§ IV–V](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevD.37.3107/fulltext)

**Decision rule.** The frozen comparison scale $|{\rm center}|+5\,{\rm uncertainty}=5.9\times10^{-21}$ is roughly $1.27\times10^{20}$ times smaller than the minimum absolute noncanonical charge. This scale is a conservative comparison, with no far-tail significance inference. A measurement incompatible with a particular exact $Q_n(k)$ rejects that branch under the stated species and calibration assumptions. A reproducible nonzero monopole incompatible with every allowed value would reject the entire declared class. The existing benchmark leaves $k=0$ empirically. [Frozen criterion](FROZEN_PROTOCOL.md) · [Exact outputs](RESULTS.json)

**Strongest comparator.** Standard gauge EFT with the observed matter representations gives $Q_n=0$, the same prediction as the selected empirical member. The qualitative no-global-phase criterion admits every branch, so it explains no previously independent measured quantity and predicts no nonzero shift against that comparator. An independently demonstrated dimension-five $B-L$ violating electroweak operator or a bare singlet Majorana mass would select $k=0$ under the frozen spectrum and interactions by gauge neutrality. The existence and operator identification of such a term would be additional physical evidence. It is unnecessary to assume either term in analyzing the known neutron result. [Novelty audit](NOVELTY_AUDIT.md)

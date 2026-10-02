# SOURCE LEDGER — GRUT SELECTOR SCREEN 01

## Grading rule and environment limitation

| grade | meaning in this screen |
|---|---|
| **PRIMARY / SOURCE-TEXT VERIFIED** | full text read. **Not achieved for any source in this screen**: arXiv / publisher / journal hosts are blocked by the egress proxy |
| **PRIMARY ABSTRACT VERIFIED (search)** | the search tool returned near-verbatim abstract text, attributed to the primary paper. Not fetched directly, so wording is close but unchecked |
| **SOURCE LOCATED, TEXT NOT RE-READ** | bibliographic record plus a partial snippet / summary only |
| **REVIEW / SECONDARY** | a review or third-party summary, used for discovery |
| **MEMORY ONLY — NOT VERDICT-BEARING** | background knowledge, not confirmed by retrieval |

**Rules applied:**
- Kill verdicts rest on **input accounting** (what the principle must be given). That accounting is visible at abstract
  level and in the frozen record.
- **No GO verdict rests on any search summary.** None was issued.
- **Retrieval route:** a census by four research sub-agents plus direct searches, all via search-engine output. arxiv.org,
  link.springer.com, quantum-journal.org, philsci-archive, OSTI and the metadata APIs were blocked.

## A. Hamiltonian / quantum mereology

| ref | source | grade | key abstract-level content used |
|---|---|---|---|
| A-1 | Carroll & Singh, "Quantum Mereology: Factorizing Hilbert Space into Subsystems with Quasi-Classical Dynamics", PRA 103, 022213 (2021), arXiv:2005.12938 | SOURCE LOCATED (snippet close to the abstract) | the factorization is chosen by a quasiclassicality criterion (an extension of the predictability sieve), given only H. A uniqueness theorem is **not** confirmed |
| A-2 | Cotler, Penington, Ranard, "Locality from the Spectrum", CMP (2019), arXiv:1702.06142 | PRIMARY ABSTRACT VERIFIED (search) | the spectrum "almost always encodes a unique description of local degrees of freedom when such a description exists"; "multiple dual local descriptions" in special cases |
| A-3 | Zanardi, "Virtual Quantum Subsystems", PRL 87, 077901 (2001) | PRIMARY ABSTRACT VERIFIED (search) | accessible observables "select a preferred tensor product structure"; compoundness "relativized" |
| A-4 | Zanardi, Lidar, Lloyd, PRL 92, 060402 (2004) | PRIMARY ABSTRACT VERIFIED (search) | TPS "relative and observable induced" |
| A-5 | Zanardi et al., "Operational Quantum Mereology and Minimal Scrambling", Quantum 8, 1406 (2024), arXiv:2212.14340 | SOURCE LOCATED | subsystems chosen by minimizing a scrambling rate over a family of algebras |
| A-6 | Cao, Carroll, Michalakis, "Space from Hilbert Space", PRD 95, 024031 (2017), arXiv:1606.08444 | SOURCE LOCATED | "Given a decomposition of Hilbert space into a tensor product of factors" |
| A-7 | Stoica, arXiv:2102.08620; arXiv:2103.15104 | SOURCE LOCATED | no-go: emergent structure that is physically relevant is not unique ("Hilbert-space fundamentalism") |
| A-8 | Soulas, Franzmann, Di Biagio, arXiv:2512.07468 (2025); Stoica comment arXiv:2603.07674 (2026) | SOURCE LOCATED | proposes a unique TPS from unitary-invariant data; contested |
| A-9 | Adil et al., "Search for classical subsystems in quantum worlds", PRD 113, 103535 (2026), arXiv:2403.10895 | SOURCE LOCATED | **several factorizations of one H admit a quasiclassical description** |
| A-10 | Loizeau & Sels, Found. Phys. 55, 3 (2025), arXiv:2409.01391 | SOURCE LOCATED | subsystem decomposition ↔ spectral decomposition, given H, an initial state and a TPS |

## B. Decoherence / Darwinism / SBS

| ref | source | grade | content used |
|---|---|---|---|
| B-1 | Zurek, Habib, Paz, PRL 70, 1187 (1993) | PRIMARY ABSTRACT VERIFIED (search) | predictability sieve selects coherent states, given the model, the split and the entropy measure |
| B-2 | Zurek, "Quantum Darwinism", Nat. Phys. 5, 181 (2009), arXiv:0903.5082 | SOURCE LOCATED | redundancy across environment fragments; threshold δ |
| B-3 | Horodecki, Korbicz, Horodecki, PRA 91, 032122 (2015) | SOURCE LOCATED | SBS given a multipartite system / environment split |
| B-4 | Riedel, PRL 118, 120402 (2017), arXiv:1608.05377 | PRIMARY ABSTRACT VERIFIED (search) | "Assuming only the tensor structure associated with spatial locality"; "A maximum length scale for records is enough to guarantee uniqueness" |
| B-5 | Kastner, arXiv:1406.4126 | PRIMARY ABSTRACT VERIFIED (search) | a distinguishable environment is "the hidden premise that makes the derivation of einselection circular" |

## C. Histories

| ref | source | grade | content used |
|---|---|---|---|
| C-1 | Dowker & Kent, PRL 75, 3038 (1995); J. Stat. Phys. 82, 1575 (1996) | PRIMARY ABSTRACT VERIFIED (search) for the PRL | "one cannot recover the standard predictions … of quasiclassical physics using the criterion of consistency alone" |
| C-2 | Gell-Mann & Hartle, PRA 76, 022104 (2007) | PRIMARY ABSTRACT VERIFIED (search) | "various coarse-grained descriptions some of which are mutually incompatible"; realms defined via chosen conserved quantities |
| C-3 | Gell-Mann & Hartle, PRA 85, 062120 (2012) | PRIMARY ABSTRACT VERIFIED (search) | a posited preferred fine-grained history basis |
| C-4 | Gell-Mann & Hartle, gr-qc/9404013; Hartle, Found. Phys. 41, 982 (2011) | SOURCE LOCATED | multiplicity of realms open |

## D. Quantum reference frames

| ref | source | grade | content used |
|---|---|---|---|
| D-1 | Giacomini, Castro-Ruiz, Brukner, Nat. Commun. 10, 494 (2019) | SOURCE LOCATED | frame-relative entanglement / superposition |
| D-2 | Ali Ahmad, Galley, Höhn, Lock, Smith, "Quantum Relativity of Subsystems", PRL 128, 170401 (2022), arXiv:2103.01232 | PRIMARY ABSTRACT VERIFIED (search) | "different reference frame perspectives induce different sets of subsystem observable algebras … frame-dependent notion of subsystems" |
| D-3 | Castro-Ruiz & Oreshkov, Commun. Phys. (2025), arXiv:2110.13199 | PRIMARY ABSTRACT VERIFIED (search) | "a quantum reference frame viewpoint is a preferred partition …; a transformation between QRFs is a change of preferred partition" |

## E–I. MaxEnt, IB, algorithmic, variational, RG

| ref | source | grade | content used |
|---|---|---|---|
| E-1 | Jaynes, Phys. Rev. 106, 620 (1957) | SOURCE LOCATED | inference given constraints |
| E-2 | Shore & Johnson, IEEE Trans. IT-26, 26 (1980) | PRIMARY ABSTRACT VERIFIED (search) | MaxEnt / min cross-entropy uniquely consistent **given constraints and a prior** |
| E-3 | Pressé, Ghosh, Lee, Dill, RMP 85, 1115 (2013) | PRIMARY ABSTRACT VERIFIED (search) | MaxCal over paths given constraints |
| E-4 | Rigol et al., PRL 98, 050405 (2007) | PRIMARY ABSTRACT VERIFIED (search) | GGE "carries more memory of the initial conditions" |
| E-5 | Ilievski et al., PRL 115, 157201 (2015) | SOURCE LOCATED + MEMORY (local-charge GGE incomplete) | constraint-set dependence of GGE (memory-grade detail, not verdict-bearing alone) |
| E-6 | Dewar 2003; Grinstein & Linsker 2007; Bruers 2007; Dewar 2009 | PRIMARY ABSTRACT VERIFIED (search) | MEP derivation contested; "not a physical principle" (Dewar 2009) |
| F-1 | Tishby, Pereira, Bialek, physics/0004057 | PRIMARY ABSTRACT VERIFIED (search) | relevance variable Y, trade-off |
| F-2 | Kolchinsky, Tracey, Van Kuyk, ICLR 2019 | PRIMARY ABSTRACT VERIFIED (search) | IB degeneracies (trivial solutions at every point when Y = f(X)) |
| F-3 | Shalizi & Crutchfield, J. Stat. Phys. 104, 817 (2001) | PRIMARY ABSTRACT VERIFIED (search) | ε-machine unique, minimal, maximally predictive **for a given process** |
| F-4 | Shalizi & Moore, "What is a macrostate?", Found. Phys. 55 (2025), cond-mat/0303625 | PRIMARY ABSTRACT VERIFIED (search) | macrostates = "the unique maximal partition" consistent with **given observables** and Markovian |
| F-5 | Koch-Janusz & Ringel, Nat. Phys. 14, 578 (2018); Gordon et al., PRL 126, 240601 (2021); Lenggenhager et al., PRX 10, 011037 (2020) | PRIMARY ABSTRACT VERIFIED (search) | IB / RSMI coarse-graining given block / buffer / environment geometry and compression level |
| G-1 | Müller, "Stationary algorithmic probability", TCS 411, 113 (2010) | PRIMARY ABSTRACT VERIFIED (search) | a machine-independent stationary prior does not exist (the attempt fails) |
| G-2 | Leike & Hutter, COLT 2015, arXiv:1510.04931 | PRIMARY ABSTRACT VERIFIED (search) | universal agent / prior is UTM-relative |
| G-3 | Wood, Sunehag, Hutter, arXiv:1111.3854 | PRIMARY ABSTRACT VERIFIED (search) | universal priors equivalent only up to constants |
| G-4 | Müller, "Law without law", Quantum 4, 301 (2020), arXiv:1712.01826 | SOURCE LOCATED | universal-machine prior over observer states |
| G-5 | Rissanen, Automatica 14, 465 (1978) | SOURCE LOCATED | description length depends on model class / code |
| H-1 | Friston, Nat. Rev. Neurosci. 11, 127 (2010); Biehl et al., Entropy 23, 293 (2021); Aguilera et al., Phys. Life Rev. 40 (2022) | PRIMARY ABSTRACT VERIFIED (search) | generative model + Markov blanket supplied; blanket conditions narrow / non-equivalent |
| I-1 | Fisher, RMP 70, 653 (1998); Latorre & Morris, JHEP 11 (2000) 004; van Enter, Fernández, Sokal, J. Stat. Phys. 72, 879 (1993) | PRIMARY ABSTRACT VERIFIED (search) | universality within a basin; scheme independence only of universal data; RG maps can be non-Gibbsian |

## J, N, O, P, Q, R(time). QES, causal sets / CDT, bootstrap, local covariance, anthropics, relational time

| ref | source | grade | content used |
|---|---|---|---|
| J-1 | Ryu–Takayanagi PRL 96, 181602 (2006); Hubeny–Rangamani–Takayanagi JHEP 0707:062; Engelhardt–Wall JHEP 01 (2015) 073 | PRIMARY ABSTRACT VERIFIED (search) | surface for a **specified boundary region** and state |
| N-1 | Rideout & Sorkin, PRD 61, 024002 (2000) | PRIMARY ABSTRACT VERIFIED (search) | a "very general family" of growth dynamics (couplings t_n) |
| N-2 | Kleitman & Rothschild, Trans. AMS 205 (1975) | SOURCE LOCATED + MEMORY (3-layer dominance) | counting dominated by non-manifold-like orders |
| N-3 | Loomis & Carlip arXiv:1709.00064; Benincasa & Dowker PRL 104, 181301 (2010) | PRIMARY ABSTRACT VERIFIED (search) | dimension-specific action; partial suppression only |
| N-4 | Ambjørn, Jurkiewicz, Loll, PRL 93, 131301 (2004); PRL 95, 171301 (2005) | PRIMARY ABSTRACT VERIFIED (search) for 2004 | 4D emergence with 4-simplices and a causal foliation as input |
| O-1 | Kos, Poland, Simmons-Duffin, JHEP 1411 (2014) 109; Kos et al., JHEP 1608 (2016) 036 | PRIMARY ABSTRACT VERIFIED (search) | island "assuming that σ and ε are the only relevant scalars" |
| O-2 | Paulos et al., arXiv:1607.06109 | PRIMARY ABSTRACT VERIFIED (search) | bounds, not a unique theory |
| P-1 | Fewster & Verch, AHP 13, 1613 (2012), arXiv:1106.4785; BFV CMP 237, 31 (2003) | PRIMARY ABSTRACT VERIFIED (search) for SPASs / dynamical locality. The **"no natural state"** theorem is **REVIEW / SECONDARY** (via arXiv:1502.04642 and the author list). BFV only "early arguments" | no locally covariant preferred state for dynamically local theories (secondary-grade) |
| P-2 | Hollands & Wald, CMP 223, 289 (2001) | PRIMARY ABSTRACT VERIFIED (search) | Wick polynomials unique up to finitely many parameters (fields, not states) |
| Q-1 | Freivogel, CQG 28, 204007 (2011); Dyson–Kleban–Susskind JHEP 0210 (2002) 011; Albrecht & Sorbo PRD 70, 063528 (2004) | PRIMARY ABSTRACT VERIFIED (search) | measure dependence; opposite conclusions from different counting |
| Q-2 | Tegmark, "On the dimensionality of spacetime", CQG 14, L69 (1997) | PRIMARY ABSTRACT VERIFIED (search) | 3+1 selected by observer existence in an ensemble ("dead worlds") |
| R-1 | Page & Wootters PRD 27, 2885 (1983); Albrecht & Iglesias PRD 77, 063506 (2008) | PRIMARY ABSTRACT VERIFIED (search) | clock and split supplied; **clock ambiguity** no-go |
| R-2 | Connes & Rovelli CQG 11, 2899 (1994) | PRIMARY ABSTRACT VERIFIED (search) | time flow from a state; outer class state-independent (frozen QFT-SCOUT G2 already priced this) |
| R-3 | Deutsch, "Constructor theory", Synthese 190 (2013) | SOURCE LOCATED | meta-framework; no selection claim |

## K, L. Quantum cosmology and arrows

| ref | source | grade | content used |
|---|---|---|---|
| K-1 | Hartle & Hawking PRD 28, 2960 (1983) | PRIMARY ABSTRACT VERIFIED (search, partial) | ground state via a postulated class of compact Euclidean geometries |
| K-2 | Hawking PRD 32, 2489 (1985); Page PRD 32, 2496 (1985); Hawking–Laflamme–Lyons PRD 47, 5342 (1993) | Page: PRIMARY ABSTRACT VERIFIED (search, summary); others SOURCE LOCATED | the CPT-symmetric total state does not fix the arrow of a WKB branch |
| K-3 | Hartle & Hertog PRD 85, 103524 (2012) | PRIMARY ABSTRACT VERIFIED (search) | two-way arrows in bouncing histories, given the no-boundary state |
| K-4 | Vilenkin PRD 30, 509 (1984); PRD 33, 3560 (1986) | 1984 PRIMARY ABSTRACT VERIFIED (search); 1986 SOURCE LOCATED | outgoing-wave boundary condition (positive frequency in a) |
| K-5 | Feldbrugge, Lehners, Turok PRD 95, 103508; PRL 119, 171301 (2017); PRD 97, 023509 (2018); Diaz Dorronsoro et al. PRD 96, 043505 (2017); Halliwell & Louko PRD 39, 2206 (1989) | PRIMARY ABSTRACT VERIFIED (search) for FLT and DDHHHJ | the lapse-contour choice (0,∞) vs (−∞,∞) changes the wavefunction |
| K-6 | Boyle, Finn, Turok PRL 121, 251301 (2018); comment arXiv:1902.07584 | PRIMARY ABSTRACT VERIFIED (search) | CPT-symmetric state imposed at the bang; continuation non-unique (comment) |
| K-7 | Carroll & Chen hep-th/0410270; Wald gr-qc/0507094; Nikolić hep-th/0411115 | PRIMARY ABSTRACT VERIFIED (search) | arrows outward on both sides; special initial conditions needed (Wald) |
| K-8 | Barbour, Koslowski, Mercati PRL 113, 181101 (2014); Zeh arXiv:1601.02790 | BKM SOURCE LOCATED; Zeh PRIMARY ABSTRACT VERIFIED (search) | Janus point with two arrows (E = L = 0); "very improbable selection condition" (Zeh) |
| K-9 | Penrose WCH (1979); Albert, *Time and Chance* (2000) | SOURCE LOCATED | arrow / low entropy postulated |
| K-10 | DeWitt, Phys. Rev. 160, 1113 (1967) | SOURCE LOCATED | boundary condition at a = 0 (exact wording unconfirmed) |
| L-1 | Bojowald PRL 87, 121301 (2001) | PRIMARY ABSTRACT VERIFIED (search) | dynamical law + initial conditions from one difference equation; **plus pre-classicality** → a unique wavefunction |
| L-2 | Cartin & Khanna PRL 94, 111302 (2005) | PRIMARY ABSTRACT VERIFIED (search) | Bianchi I: "only the zero solution" satisfies the constraints |
| L-3 | Date PRD 72, 067301 (2005); Bojowald, Cartin, Khanna PRD 76, 064018 (2007); Bojowald & Simpson CQG 31, 185016 (2014) | BKS / BS PRIMARY ABSTRACT VERIFIED (search); Date SOURCE LOCATED | lattice-refinement and factor-ordering dependence |

## R. Other candidates found

| ref | source | grade | content used |
|---|---|---|---|
| R-4 | Brandenberger & Vafa, NPB 316, 391 (1989) | SOURCE LOCATED | winding modes annihilate generically only in ≤ 3 large spatial dimensions |
| R-5 | Greene, Kabat, Marnerides, PRD 88, 043527 (2013), arXiv:1212.2115 | PRIMARY ABSTRACT VERIFIED (search) | "three is the maximum number … from an initial thermal fluctuation"; preferred three "if the string coupling is sufficiently large"; "a more careful study … is needed to assess the likelihood of the assumptions" |
| R-6 | Easther, Greene, Jackson, Kabat, JCAP 0502 (2005) 009, hep-th/0409121 | SOURCE LOCATED | string-windings follow-up (content not retrieved) |

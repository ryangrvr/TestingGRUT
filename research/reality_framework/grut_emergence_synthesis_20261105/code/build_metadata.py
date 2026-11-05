"""Build the human/machine claim ledger, bibliography and a schematic map."""
import json
import textwrap
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = [
 ("S01","Donoghue (1994), General relativity as an effective field theory",
  "https://arxiv.org/abs/gr-qc/9405057","Primary abstract; low-energy scope"),
 ("S02","Donnelly and Freidel (2016), Local subsystems in gauge theory and gravity",
  "https://arxiv.org/abs/1601.04744","Primary abstract; boundary/subsystem qualification"),
 ("S03","Page and Wootters (1983), Evolution without evolution",
  "https://doi.org/10.1103/PhysRevD.27.2885","Publisher abstract; relational-clock proposal"),
 ("S04","Schlosshauer (2019), Quantum Decoherence",
  "https://arxiv.org/abs/1911.06282","Primary review abstract; not a new measurement rule"),
 ("S05","Reeb and Wolf (2014), An improved Landauer Principle with finite-size corrections",
  "https://arxiv.org/abs/1306.4352","Primary abstract; finite equality rederived in this packet"),
 ("S06","D'Onofrio and Rummukainen (2016), The Standard Model cross-over on the lattice",
  "https://arxiv.org/abs/1508.07161","Primary abstract and reported crossover temperature"),
 ("S07","CERN, Heavy ions and quark-gluon plasma",
  "https://home.cern/science/physics/heavy-ions-and-quark-gluon-plasma","Institutional account of early hadronic matter"),
 ("S08","NASA Science, Overview of the universe",
  "https://science.nasa.gov/universe/overview/","Nucleosynthesis/recombination passages; not adopting all inflation wording"),
 ("S09","NASA Science, What were the first stars like?",
  "https://science.nasa.gov/mission/webb/science-overview/science-explainers/what-were-the-first-stars-like/",
  "First-star timing passage; uncertain chronology"),
 ("S10","Kolchinsky and Wolpert (2018), Semantic information, autonomous agency and non-equilibrium statistical physics",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6227811/","Intervention/viability sections; theoretical framework, not consciousness criterion"),
 ("S11","COGITATE Consortium et al. (2025), Adversarial testing of global neuronal workspace and integrated information theories of consciousness",
  "https://doi.org/10.1038/s41586-025-08888-1","Primary indexed abstract/discussion; direct publisher opening encountered a redirect"),
 ("S12","Burgess (2013), The Cosmological Constant Problem",
  "https://arxiv.org/abs/1309.4133","Primary review abstract; EFT naturalness scope"),
 ("S13","Klinkhamer and Volovik (2008), Self-tuning vacuum variable and cosmological constant",
  "https://arxiv.org/abs/0711.3170","Primary abstract; imperfect-vacuum interpretation credited"),
 ("S14","Klinkhamer and Volovik (2008), Dynamic vacuum variable and equilibrium approach in cosmology",
  "https://arxiv.org/abs/0806.2805","Primary abstract; a specific proposed dynamics"),
 ("S15","Chandrasekaran, Longo, Penington and Witten (2023), An Algebra of Observables for de Sitter Space",
  "https://arxiv.org/abs/2206.10780","Primary abstract; observer-worldline dressing, not consciousness"),
 ("S16","Gibbons and Hawking (1977), Cosmological event horizons, thermodynamics, and particle creation",
  "https://doi.org/10.1103/PhysRevD.15.2738","Publisher abstract; de Sitter entropy identity algebra checked here"),
]

# Each challenge is scoped: a consistency test is not an empirical falsifier.
ROWS = [
 ("G01","ADOPTED","Quantum/Born structure and the selected dynamics are premises.",
  "Specified Hilbert/algebraic model and operational effects.","R01 R04","A weaker derivation would upgrade this; closure terminology alone does not."),
 ("G02","CONDITIONAL","No additional GRUT law was obtained in the inspected closure campaigns.",
  "Their stated model classes and searches only.","R04 R06","Exhibit an additional consequence from the stated premises; no universal no-go is claimed."),
 ("G03","CONDITIONAL","A whole and individually measurable parts can coexist.",
  "Justified tensor or algebraic subsystem description.","S02 W01","Show inconsistency of the chosen restriction or reproduce the finite countercheck."),
 ("G04","CONDITIONAL","Equal local marginals can have different local futures.",
  "Two correlated two-qubit preparations and H=XX.","W01","Recompute the trace-distance result 1/sqrt(2) at pi/8."),
 ("G05","OPEN","A unique physical subsystem/access structure is not derived here.",
  "Continuum/gauge setting must be specified.","R04 S02","Construct and validate a selection rule under weaker premises."),
 ("G06","ADOPTED","The tuple is a theory specification, not a unique solution.",
  "Its six independently declared slots.","R04","Supply an actual constraint algebra, state and observable dictionary."),
 ("G07","CONDITIONAL","A constrained stationary finite history can write and erase records.",
  "Eight-reading clock, reversible gates and a supplied blank register.","W01 S03","Recompute the propagation kernel and conditional correlations."),
 ("G08","CONDITIONAL","Maximum mixing in that finite physical subspace gives no conditional record history.",
  "The specific eight-dimensional history subspace.","W01","Compute its clock-conditioned maximally mixed state."),
 ("G09","CONDITIONAL","The Gibbs-reservoir entropy identity holds.",
  "Finite systems, initial product state, Gibbs complement, joint unitary.","S05 W01","Check the entropy, mutual information, relative entropy and total-energy residuals."),
 ("G10","CONDITIONAL","Local entropy need not increase monotonically.",
  "Finite closed example and specified subsystem.","W01","Recompute the zero, ln(2), zero cycle; not a claim against macroscopic thermodynamics."),
 ("G11","STANDARD","The Standard Model electroweak event is a smooth crossover near 160 GeV.",
  "Physical Higgs regime and Standard Model thermal history.","S06","Compare to the source's lattice result within its assumptions."),
 ("G12","STANDARD","BBN, recombination and stars occupy different, overlapping stages of cosmic history.",
  "Hot cosmology and astrophysical uncertainties.","S07 S08 S09","Compare chronology and physical regimes; no precise life/experience onset supplied."),
 ("G13","CONDITIONAL","The zero-density M1 tree potential has one real stationary point.",
  "Listed epsilon/lambda ratio; X=0; tree level.","R07","Factor its derivative and check the quadratic discriminant -0.16."),
 ("G14","CONDITIONAL","Object charge fraction eta can exceed one with compensating negative gas.",
  "Positive cores plus unit-charge free X and anti-X; common denominator.","R07","Check global charge conservation; add missing charged species if present."),
 ("G15","CONDITIONAL","Core conversion is not automatically irreversible cleanup.",
  "Energy and charge conserving forward/reverse channels, appropriate spectra.","R07","Account for reverse thresholds and neutral decays in a complete kinetic calculation."),
 ("G16","FAILED-IN-SCOPE","Zero of 135 keys passes strict M2-X.",
  "Frozen scan rule, two simulated nodes, archived numerical closure.","R08 R09","Reaggregate saved roots with the same rule; controlled new physics/history is a new test."),
 ("G17","CONDITIONAL","Two isolated Q roots and 15 extrapolated free-X diagnostics do not repair G16.",
  "Different acceptance conditions; extrapolations excluded in advance.","R08 R09","Check bracket flags and branch-specific rules."),
 ("G18","CONDITIONAL","Fresh H55 scalar evolution changes the hand-over radiation fraction.",
  "Same transferred H84 fermion drain; one scalar resolution.","R08 R10","Compare saved controls; do not infer a new relic without a new splice."),
 ("G19","OPEN","Final actual object abundance and velocity distribution are uncertified.",
  "Full history, population extraction, controlled fermions and reactions.","R08 R10 R13","Complete connected evolution with conservation, matching and convergence."),
 ("G20","FAILED-IN-SCOPE","The four registered co-genesis candidates fail their condition.",
  "Specified source maps, textures and histories.","R11","Recompute those maps; no exclusion of all common-origin mechanisms."),
 ("G21","CONDITIONAL","An invariant unitary charge-flipping symmetry gives zero mean charge.",
  "Unitary U commutes with H and preserves rho0; U†QU=-Q.","R04","Verify the manuscript's algebraic proof; antiunitary time reversal is a different statement."),
 ("G22","OPEN","Cosmic state/phase selection is not supplied by closure.",
  "Within the inspected preparations and framework.","R04 R06 R12","Provide a physical selection rule with a new consequence."),
 ("G23","CONDITIONAL","Information scrambling lowers a supplied controller's success from .9 to .5.",
  "Specified sensor, environment, policy and success function.","W01 S10","Repeat the matched intervention; controller origin and energetic cost remain separate."),
 ("G24","OPEN","Consciousness does not follow from these finite record or control examples.",
  "No bridge criterion is assumed.","S11 S15","State a discriminating bridge and test it across interventions."),
 ("G25","PROPOSED","Underlying consciousness or intention requires a separately specified hypothesis.",
  "Alternatives and observations must be defined.","S11","Give likelihoods that differ from natural physical alternatives; anomaly alone is insufficient."),
 ("G26","CONDITIONAL","A constant H shift is invisible to the fixed-background density evolution.",
  "Constant times identity; no dynamical gravitational coupling.","R05 W01","Recompute unitary density evolution; a gravity coupling leaves the scope."),
 ("G27","STANDARD","A vacuum constant gravitates as stress-energy proportional to the metric.",
  "Metric variation of -C integral sqrt(-g).","R05 S12","Check metric variation and renormalized matching convention."),
 ("G28","ADOPTED","The card's 2.52e-47 GeV^4 is a measured vacuum calibration.",
  "Archived M2 cosmological calibration, not a new measurement.","R07","Derive its sign/scale stably rather than fit the same datum."),
 ("G29","CONDITIONAL","The de Sitter entropy/vacuum formula is an identity within its assumptions.",
  "Four-dimensional Einstein de Sitter gravity and reduced Planck units.","S16","Derive rho=24 pi^2 M_P^4/S; specify S independently to predict rho."),
 ("G30","OPEN","No GRUT-specific universal coupling relation or new distinctive signal is established here.",
  "This synthesis and inspected campaigns only.","R04 R06 R07 R08","Produce a relation removing a free input plus a frozen held-out observable."),
 ("G31","PROPOSED","A complete working account can coexist with incomplete fundamental explanation.",
  "Definition of three different finish lines.","R01 R02 R03","Audit whether any model success is being used to discharge an unrelated open problem."),
]


def main():
    m=json.loads((ROOT/"SOURCE_MANIFEST.json").read_text())
    m["external"]=[dict(id=i,title=t,url=u,read_scope=r) for i,t,u,r in EXTERNAL]
    (ROOT/"SOURCE_MANIFEST.json").write_text(json.dumps(m,indent=2)+"\n")
    out=["# Sources and provenance","",
         "Pinned project snapshot: "+m["base_commit"]+". Relevant sections were read; linked campaigns were not fully rerun.",
         "External entries distinguish primary abstracts, selected passages and full-work claims. No citation is an independent replication.","",
         "## Project evidence",""]
    for r in m["local"]:
        out += [f"- **{r['id']}** [{r['path']}]({r['url']}). Exact Git blob: {r['git_blob']}; SHA-256 in the JSON manifest."]
    out += ["","## Literature and institutional chronology",""]
    for i,t,u,r in EXTERNAL: out += [f"- **{i}** [{t}]({u}). Read scope: {r}."]
    out += ["","## Computations in this packet","",
            "- **W01** [Finite witnesses](code/witnesses.py) and [raw results](results/witnesses.json). Same-author reconstruction; known analytic answers. No empirical data or cosmological solver execution.",""]
    (ROOT/"SOURCES.md").write_text("\n".join(out))
    dependencies={"G04":["G03"],"G08":["G07"],"G17":["G16"],"G19":["G16","G18"],
                  "G22":["G02"],"G24":["G07","G23"],"G30":["G16","G22"],"G31":["G02","G19","G24","G30"]}
    claims=[]
    for i,status,text,premises,source,challenge in ROWS:
        claims.append(dict(id=i,status=status,claim=text,premises=premises,
                           sources=source.split(),challenge=challenge,depends_on=dependencies.get(i,[]),
                           independent_review_in_this_packet=False))
    (ROOT/"CLAIMS.json").write_text(json.dumps(claims,indent=2)+"\n")
    out=["# Claim ledger","", "Statuses apply to the stated scope. New calculations have no second-author review.",
         "The challenge column contains different kinds of test: mathematical, numerical, conceptual or empirical. They are not interchangeable.","",
         "| ID / status | Claim | Premises | Challenge | Sources |","|---|---|---|---|---|"]
    for c in claims:
        out += [f"| {c['id']} / {c['status']} | {c['claim']} | {c['premises']} | {c['challenge']} | {', '.join(c['sources'])} |"]
    (ROOT/"CLAIMS.md").write_text("\n".join(out)+"\n")
    make_map()
    print(json.dumps({"project_sources":len(m["local"]),"external_sources":len(EXTERNAL),"claims":len(claims)}))


def make_map():
    # A dependency diagram, deliberately not a numeric cosmic-time plot.
    nodes=[
      ("laws",.25,.90,"Quantum and gravity framework","ADOPTED","R04 S01",["Fields, symmetries, observation rule"]),
      ("state",.75,.90,"Initial state and orientation","OPEN","R04",["A preparation is supplied"]),
      ("early",.50,.76,"Expansion and early interactions","CONDITIONAL","R07",["Inflation/reheating: model dependent"]),
      ("ordinary",.25,.60,"Ordinary matter history","STANDARD","S06 S07 S08",["Electroweak / QCD / nuclei / atoms"]),
      ("dark",.75,.60,"Candidate dark sector","OPEN","R08 R09",["Strict M2-X: 0 of 135 passes"]),
      ("structure",.50,.44,"Persistent physical structures","STANDARD","S09",["Gravity, binding, chemistry, energy flows"]),
      ("records",.25,.28,"Records and correlations","CONDITIONAL","W01 S04",["Finite mechanism demonstrated"]),
      ("agency",.75,.28,"Information use and feedback","CONDITIONAL","W01 S10",["Prepared controller; origin not derived"]),
      ("experience",.50,.11,"Subjective experience","OPEN","S11",["A physical bridge remains missing"]),
    ]
    edges=[("laws","early","dynamics",False),("state","early","preparation",True),
           ("early","ordinary","thermal history",False),("early","dark","candidate history",True),
           ("ordinary","structure","assembly",False),("dark","structure","possible contribution",True),
           ("structure","records","interactions",False),("structure","agency","organization",True),
           ("records","agency","usable information",False),("agency","experience","unknown bridge",True)]
    payload={"meaning":"Dependency graph, not an evenly spaced timeline. Branches overlap; open arrows are not demonstrated.",
             "cross_cutting_open":["origin of laws/state","vacuum density","experience"],
             "nodes":[dict(id=i,x=x,y=y,label=l,status=s,sources=refs.split()) for i,x,y,l,s,refs,_ in nodes],
             "edges":[dict(source=a,target=b,label=l,conditional_or_open=d) for a,b,l,d in edges]}
    (ROOT/"EMERGENCE_MAP.json").write_text(json.dumps(payload,indent=2)+"\n")
    fig,ax=plt.subplots(figsize=(10,12),dpi=160)
    bg="#f7f7f1"; fig.patch.set_facecolor(bg); ax.set_facecolor(bg)
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    fig.subplots_adjust(left=.03,right=.97,top=.92,bottom=.055)
    fig.text(.055,.965,"GRUT · emergence and interaction",fontsize=22,fontweight="bold",color="#173e49")
    fig.text(.055,.935,"One physical whole; interacting parts with distinct histories",fontsize=12,color="#425760")
    colors={"ADOPTED":"#e7e6ed","OPEN":"#f8e6d1","CONDITIONAL":"#e2edf5","STANDARD":"#e1eee4"}
    positions={n[0]:(n[1],n[2]) for n in nodes}; width=.39;height=.085
    for a,b,label,dashed in edges:
        xa,ya=positions[a];xb,yb=positions[b]
        if abs(ya-yb)<.01:
            start=(xa+width/2+.009,ya);end=(xb-width/2-.013,yb)
        else:
            start=(xa,ya-height/2-.009);end=(xb,yb+height/2+.014)
        ax.add_patch(FancyArrowPatch(start,end,arrowstyle="-|>",mutation_scale=12,
                                    linewidth=1.2,linestyle="--" if dashed else "-",
                                    color="#6c7d83",zorder=1))
        if a=="records" and b=="agency":
            ax.text(.5,.295,"information",ha="center",va="bottom",fontsize=8,color="#566b75")
    for i,x,y,l,s,refs,detail in nodes:
        ax.add_patch(FancyBboxPatch((x-width/2,y-height/2),width,height,
                                  boxstyle="round,pad=0.008,rounding_size=0.008",
                                  edgecolor="#92a4a7",facecolor=colors[s],linewidth=.8,zorder=3))
        ax.text(x,y+.018,l,fontsize=11.0,weight="bold",ha="center",va="center",color="#173e49",zorder=4)
        ax.text(x,y-.006,detail[0],fontsize=9,ha="center",va="center",color="#425760",zorder=4)
        ax.text(x,y-.028,s+"  ·  "+refs.replace(" ",", "),fontsize=8,ha="center",va="center",color="#566b75",zorder=4)
    fig.text(.055,.037,"Schematic dependencies, not a scaled chronology. Dashed paths require model choices or an open bridge.",
             fontsize=9,color="#425760")
    fig.text(.055,.020,"Vacuum density and the origin of laws remain open across the map. Source snapshot: 8ab5829a.",
             fontsize=9,color="#425760")
    fig.savefig(ROOT/"EMERGENCE_MAP.svg",facecolor=bg)
    fig.savefig(ROOT/"EMERGENCE_MAP.png",facecolor=bg)
    plt.close(fig)


if __name__ == "__main__": main()

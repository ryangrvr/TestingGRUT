"""Assemble the theory datum, verification summary and source manifest."""
import ast
import hashlib
import json
import re
import subprocess
from common import ROOT,write,read_previous

REPO=ROOT.parents[2]
PARENT='9877518b65989c79adb742e1ff33d7d4b0f85a7d'


def main():
    a=json.loads((ROOT/'HISTORY_RESULTS.json').read_text())
    b=json.loads((ROOT/'OPERATOR_RESULTS.json').read_text())
    groups={**a['groups'],**b['groups']}
    assert len(groups)==8 and all(g['passed'] for g in groups.values())
    write('VERIFICATION_RESULTS.json',dict(groups=groups,passed=8,total=8,
          physical_cosmology_certified=False,quantum_gravity_derived=False,closure_state_law_derived=False))
    previous=read_previous()
    write('THEORY_DATA.json',dict(
        architecture='Constrained relational quantum description with covariant EFT; no fundamental completion selected.',
        spacetime=dict(signature='-+++',spatial_branch='T3, supplied compact topology',gravity='Adopted Einstein EFT plus allowed counterterms'),
        fields=dict(Standard_Model='Adopted',binding_phi='real, hidden neutral',X=1,Y=3,Sigma=6,gauge_group='SM x U(1)_X',reporter_D='Optional earlier probe, not required by this core architecture'),
        scalar_basis=dict(max_mass_dimension=4,monomials=26,with_exact_phi_parity=19,reference_potential='THEORY.md; full renormalizable basis ACTION.md'),
        matching=dict(low_energy_targets=previous['light_sector'],heavy_reference=previous['normalized_parameters'],
                      missing='Complete quantum matching, higher-dimensional coefficients and physical inflation/phase-holder state.',rho_vac='Independent renormalized parameter'),
        additional_curvature_representative=dict(xi_X=-5/6,mass_squared_over_H2=dict(de_Sitter=-10,matter=-2.5,radiation=0),
             status='Added covariant representative; no new cosmological history computed.'),
        observations=dict(probability='Tr_phys(rho_phys E_a_given_t)',Born_rule='Supplied',finite_dictionary='RECORDS.md',continuum_gravity_dictionary='Unconstructed'),
        exact_finite_instance=dict(kinematic_dimension=64,physical_dimension=8,clock_readings=8,input_state='|000> at clock reading0',
                                  phase=.37,record_window=[3,4,5],state_preparation='Supplied'),
        predictions=dict(kinetic_solver='Claude, existing workflow',experimental_family='Previously selected multi-target nuclear recoil shapes',distinctive_signal_established=False),
        missing_fundamental_selection=dict(action_or_state_rule='No new restriction derived from closure',vacuum_density=False,cosmic_arrow=False,subjective_consciousness=False)))

    sources=[]
    specs=[
      ('grut_closed_whole_theory_20261019/FORMAL_CORE.md','Definitions and theory datum, approximately lines1–195; state-selection discussion lines415–495. Not executed.'),
      ('grut_cstar_audit_20261024/CSTAR_HYPOTHESIS.md','Tagged hypothesis, distinctions between state spaces, constraints, topology and vacuum. Read, not rerun.'),
      ('grut_concept_checkpoint_20261025/FOUNDATION_GATE.md','Full principle-level restriction requirement and unfinished foundations. Read.'),
      ('grut_concept_checkpoint_20261025/REPORT.md','Concept status, empirical scope and source/debris distinction. Read.'),
      ('grut_z6_gauge_bridge_20261025/CONSTRUCTION.md','Full tree action/matching declaration and explicit open quantum portals/hierarchy. Read.'),
      ('grut_heavy_vector_budget_20261025/INPUTS.json','Numerical supplied point and preparation metadata. Read by operators.py; no release or mode integration rerun.'),
      ('grut_kinetic_contract_20261026/KINETIC_SPEC.md','Scope, state/formation and quantum-kernel assumptions. Read; solver remains with Claude.')
    ]
    for suffix,scope in specs:
        path='research/reality_framework/'+suffix
        d=subprocess.check_output(['git','show',PARENT+':'+path],cwd=REPO)
        sources.append(dict(kind='git_file',commit=PARENT,path=path,bytes=len(d),sha256=hashlib.sha256(d).hexdigest(),read_scope=scope))
    sources.extend([
      dict(kind='primary_paper',url='https://journals.aps.org/prd/abstract/10.1103/PhysRevD.27.2885',read_scope='Official abstract only. Relational-time prior art; no full-text access claim.'),
      dict(kind='primary_paper',url='https://arxiv.org/pdf/quant-ph/0405098',read_scope='History-state construction page4 and input/propagation Hamiltonians pages9–10. Cyclic variant and checks are ours.'),
      dict(kind='primary_paper',url='https://arxiv.org/html/hep-ph/0111209v2',read_scope='Introduction loop normalization and section3, eqs3.2–3.4. One-loop Landau/MS-bar potential; no two-loop result evaluated.')
    ])
    write('SOURCES.json',dict(records=sources,new_constructions='Exact finite history/POVM example; scalar rational polynomial; scoped gauge derivative; curvature representative.'))
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or p.name=='MANIFEST.json' or '__pycache__' in p.parts:continue
        d=p.read_bytes()
        if p.suffix=='.py':ast.parse(d,filename=str(p))
        if p.suffix=='.json':json.loads(d,parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
        if p.suffix=='.md':
            for dest in re.findall(r'\[[^\]]*\]\(([^)]+)\)',d.decode()):
                if dest.startswith(('http://','https://','#')):continue
                assert (p.parent/dest.split('#')[0]).exists(),(p,dest)
        files.append(dict(path=str(p.relative_to(ROOT)),bytes=len(d),sha256=hashlib.sha256(d).hexdigest()))
    write('MANIFEST.json',dict(files=files,excludes=['MANIFEST.json','__pycache__'],source_records=len(sources)))
    print('PACKET_OK',len(files),'files',len(sources),'source records; mathematical checks8/8')


if __name__=='__main__':main()

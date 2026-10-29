"""Provenance, conservation, and control comparisons for 3D multipoles."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DYN=ROOT.parent/'grut_m1_dynamics_20261027'

def main():
    parent=json.loads((DYN/'RUNS.json').read_text())['cases'];archive=np.load(DYN/'TRACE_DATA.npz')
    rows=[]
    for name,oldname in [('object64_control','object64'),('object64_perturbed','object64_perturbed'),
                         ('object96_control',None),('object96_perturbed','object96_perturbed'),
                         ('object96_control_double',None),('object96_perturbed_double','object96_perturbed')]:
        path=ROOT/'.work'/name;trace=np.genfromtxt(path/'state.tsv',names=True,delimiter='\t')
        initial=hashlib.sha256((path/'initial.bin').read_bytes()).hexdigest()
        row={'name':name,'initial_sha256':initial,
             'max_relative_energy_drift':float(np.max(abs(trace['rho']/trace['rho'][0]-1))),
             'max_relative_charge_drift':float(np.max(abs(trace['nQ']/trace['nQ'][0]-1))),
             'last_time':float(trace['t'][-1])}
        if oldname:
            old=parent[oldname];row['matches_archived_initial_sha256']=initial==old['parameters']['initial_sha256']
            ref=archive[oldname];comparison=[]
            for t in (0.,26.,52.,78.):
                if t>trace['t'][-1]:continue
                a=trace[np.argmin(abs(trace['t']-t))];b=ref[np.argmin(abs(ref['t']-t))]
                if abs(a['t']-b['t'])>1e-9:raise ValueError('No common trace time')
                comparison.append({'t':t,'relative_energy_difference':float(a['rho']/b['rho']-1),
                                   'relative_charge_difference':float(a['nQ']/b['nQ']-1)})
            row['common_time_comparison_to_archived_trace']=comparison
        rows.append(row)
    studies={}
    for name in ('multipoles64','multipoles96','multipoles96_ang12','multipoles96_dr05',
                 'multipoles96_l6','multipoles64_l6','multipoles96_t20','multipoles96_double',
                 'multipoles96_double_l6','multipoles96_double_dr05','multipoles96_double_t20'):
        d=json.loads((ROOT/(name+'.json')).read_text())['result']
        studies[name]={'last_time':d['last_time'],'orthonormality_max_error':d['orthonormality_max_error'],
                       'r20_difference':[r for r in d['rows'] if r['radius_mX_inverse']==20 and r['case']=='difference_to_control']}
    out={'3D_runs':rows,'multipole_controls':studies,
         'scope':'Periodic side-64 box, 5% disturbance, early finite-radius field multipoles. Archived 96^3 full field snapshots were not present in the parent packet; these are same-code reproductions of its initial conditions and early trace, not independent replication.'}
    (ROOT/'ANGULAR_AUDIT.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()

"""Mask-expansion sensitivity; not a background-subtracted soliton mass."""
from analyze import *

def expand(mask,steps):
    out=mask.copy()
    for _ in range(steps):
        nxt=out.copy()
        for axis in range(3):nxt|=np.roll(out,1,axis)|np.roll(out,-1,axis)
        out=nxt
    return out

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];analysis=json.loads((ROOT/'ANALYSIS.json').read_text())['cases'];out={}
    for tag in ['mixed_N96','onesign_N96','mixed_N128']:
        if tag not in analysis:continue
        row,q,p,j,labs=analyze_snapshot(runs[tag]['snapshots'][-1],details=True)
        rho=local_energy(q,p,row['L']);dv=(row['L']/row['N'])**3;rows=[]
        candidates=[c for c in row['catalogs'][1]['candidates'] if c['compact']][:3]
        for c in candidates:
            m=labs==c['id'];values=[]
            for layers in [0,1,2,4]:
                mask=expand(m,layers);Q=float(j[mask].sum()*dv);E=float(rho[mask].sum()*dv)
                values.append(dict(layers=layers,nominal_extension_m_inverse=layers*row['L']/row['N'],cells=int(mask.sum()),Q_mask=Q,E_mask_m=E,E_mask_over_Qm=E/Q))
            rows.append(dict(original_candidate=c,expansions=values))
        out[tag]=rows;save('MASK_SENSITIVITY.json',dict(cases=out,
            scope='Periodic six-neighbour dilation, not spherical radii. Regions may include ambient gas or other structures; energies have no background subtraction and are not isolated-object masses.'))
        print('EXPANDED',tag,[[v['E_mask_over_Qm'] for v in r['expansions']] for r in rows],flush=True)

if __name__=='__main__':main()

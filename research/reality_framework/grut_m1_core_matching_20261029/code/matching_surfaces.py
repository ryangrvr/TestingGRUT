"""Same-run finite-radius flux and full shell-field energy/charge ledger."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'grut_m1_local_cloud_20261029'
G2=.01

def filehash(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def shell(path,lo,hi,h=.05):
    d=np.genfromtxt(path,names=True,delimiter=',');r=d['rho'];n=len(r)
    edges=np.arange(n+1)*h;vol=4*np.pi/3*(edges[1:]**3-edges[:-1]**3)
    face=4*np.pi*edges[1:]**2/h;face[-1]*=2
    x=d['x'];y=d['y'];s=d['s'];px=d['px'];py=d['py'];ps=d['ps'];q=x*x+y*y
    diff=lambda a,b:np.diff(np.r_[a,b])
    e=vol*(.5*(px*px+py*py+ps*ps)+.75*(s*s-1)**2-.2*(s-1)**3+.5*s*s*q+.025*q*q)
    e+=.5*face*(diff(x,0)**2+diff(y,0)**2+diff(s,1)**2)
    charge=vol*(y*px-x*py);mask=(r>=lo)&(r<hi)
    return {'energy_mX':float(e[mask].sum()/G2),'charge':float(charge[mask].sum()/G2),
            'num_cells':int(mask.sum()),'r_inner':lo,'r_outer':hi}

def result_at(radius,width):
    if width==400:
        z=np.load(OLD/f'source_h005_r{radius}.npz');d=json.loads((OLD/f'source_h005_r{radius}.json').read_text())
    else:
        z=np.load(ROOT/'.work'/f'probe_r{radius}_t{width}.npz');d=json.loads((ROOT/'.work'/f'probe_r{radius}_t{width}.json').read_text())
    row={'r':radius,'window':width,'probe_actual':d['probe_actual_mX_inverse'],
         'raw_charge':d['raw_flux']['charge'],'raw_energy_mX':d['raw_flux']['charged_energy_mX']+d['raw_flux']['neutral_energy_mX']}
    row['species']={}
    for s in ('plus','minus','neutral'):
        v=d['species'][s];pos=z[s+'_positive'];signed=z[s+'_signed']
        row['species'][s]={'positive_projection_number':v['positive_projection_number'],'signed_number':v['signed_number'],
                            'negative_bin_weight':v['negative_bin_weight'],
                            'first_energy_bin_positive_number':float(pos[:,0].sum()),
                            'first_energy_bin_signed_number':float(signed[:,0].sum())}
    return row

def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,default=OLD/'.work')
    p.add_argument('--initial',type=Path,default=ROOT.parent/'grut_m1_expansion_clock_20261028/.work/ringdown_h005_initial.csv')
    args=p.parse_args();raw=args.raw;initial=args.initial
    before=shell(initial,40,80);after=shell(raw/'radial_r40_final.csv',40,80)
    a=np.genfromtxt(raw/'radial_r40_history.csv',delimiter=',',names=True)
    b=np.genfromtxt(raw/'radial_r80_history.csv',delimiter=',',names=True)
    e_inout=float((a['Eflux_probe'][-1]-b['Eflux_probe'][-1])/G2)
    q_inout=float((a['Qflux_probe'][-1]-b['Qflux_probe'][-1])/G2)
    ledger={'initial':before,'final':after,'delta_energy_mX':after['energy_mX']-before['energy_mX'],
            'delta_charge':after['charge']-before['charge'],
            'in_minus_out_integrated_probe_energy_mX':e_inout,'in_minus_out_integrated_probe_charge':q_inout,
            'energy_residual_boundary_matching_mX':after['energy_mX']-before['energy_mX']-e_inout,
            'charge_residual_boundary_matching':after['charge']-before['charge']-q_inout,
            'interpretation':'Probe flux uses central derivatives at cell centers; shell energy uses discrete finite-volume cell and face assignment. Small residual is a discretization/matching diagnostic, not unseen new physics.'}
    out={'initial_sha256':filehash(initial),'final_sha256':filehash(raw/'radial_r40_final.csv'),
         'energy_bin_edges_first_three':[1,1.025,1.05],
         'rows':[result_at(r,w) for w in (200,400,800) for r in (40,60,80)],'field_shell_40_to_80':ledger,
         'scope':'Same h=.05 radial run, surface extraction at fixed coordinate time. These are flux and shell-storage diagnostics; no bound follows for interior interactions, waves not yet crossed by the ending time, quantum particles, or final mass.'}
    (ROOT/'MATCHING_SURFACES.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()

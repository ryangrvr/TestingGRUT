#!/usr/bin/env python3
"""Report actual grid/time/boundary evidence, without forcing extrapolation."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scattering import ROOT


def main():
    rows=[json.loads(p.read_text()) for p in sorted(ROOT.glob('[HN]*_RESULTS.json'))]
    byname={r['name']:r for r in rows}
    snapshots=[]
    for row in rows:
        name=row['name'];tab=np.loadtxt(ROOT/'data'/f'{name}_final.csv',delimiter=',',skiprows=1)
        r=tab[:,0];q=tab[:,2]*tab[:,4]-tab[:,1]*tab[:,5]
        mask=r<=40.;cum=cumulative_trapezoid(4*np.pi*r[mask]**2*q[mask],r[mask],initial=0.)
        r90=float(np.interp(.9*cum[-1],np.maximum.accumulate(cum),r[mask]))
        history=np.loadtxt(ROOT/'data'/f'{name}_history.csv',delimiter=',',skiprows=1)
        at_times={str(t):float(np.interp(t,history[:,0],history[:,12])/history[0,2])
                  for t in [1000.,2000.,4000.,8000.,12000.] if t<=row['parameters']['tmax']}
        snapshots.append(dict(name=name,core_charge_radius90_within40=r90,
            points_across_radius90=r90/row['parameters']['h'],core60_at_times=at_times))
    fine=[byname[n]['core60_late_mean_fraction'] for n in ['N025','N0125','N00625']]
    d1,d2=fine[0]-fine[1],fine[1]-fine[2]
    monotone=d1*d2>0
    order=float(np.log2(abs(d1/d2))) if monotone else None
    fixed4=[(16*fine[i+1]-fine[i])/15 for i in [0,1]]
    spatial=dict(names=['N025','N0125','N00625'],values=fine,
        relative_spread=(max(fine)-min(fine))/np.mean(fine),
        observed_errors_monotone=monotone,observed_Richardson_order=order,
        formal_fixed_order4_extrapolants=fixed4,
        extrapolation_certified=bool(monotone and order is not None and 3<order<5),
        extrapolation_caution='Fixed-order extrapolants are diagnostics; nonmonotone differences do not establish the asymptotic order.')
    controls=[]
    for a,b in [('N025','N025_time'),('N0125','N0125_time'),('N0125','N0125_abs'),
                ('N0125','N0125_box'),('N0125_box','N0125_far')]:
        if a not in byname or b not in byname:continue
        x=byname[a];y=byname[b]
        controls.append(dict(a=a,b=b,retained_charge_relative_change=y['core60_late_mean_fraction']/x['core60_late_mean_fraction']-1,
            energy_error_ratio=x['max_energy_error']/y['max_energy_error']))
    out=dict(scope='One supplied Gaussian, spherical finite-time core60 observable; no actual cosmological dN/dQ.',
        runs=rows,snapshot_diagnostics=snapshots,spatial=spatial,controls=controls,
        finite_time_one_percent_agreement=bool(spatial['relative_spread']<.01 and
            all(abs(c['retained_charge_relative_change'])<.01 for c in controls)),
        asymptotic_object_or_3D_population_certified=False)
    (ROOT/'CONVERGENCE_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(spatial=spatial,controls=controls),indent=2))


if __name__=='__main__':main()

from pathlib import Path
import json,numpy as np
p=Path(__file__).resolve();f=p.parents[1]/'cp09_bounded_clock/check.py';ns={'__file__':str(f)};src=f.read_text().split('rows=[]')[0]
# Preserve frozen simulator; declare a new pulse-calibration control only.
src=src.replace('S=expm(1j*sign*eps/width*Bk*dt);D=expm(1j*zeta/width*Ak*dt)','ws=source.sum()*dx;wd=det.sum()*dx\n S=expm(1j*sign*eps/ws*Bk*dt);D=expm(1j*zeta/wd*Ak*dt)')
src=src.replace('sign*eps/width*np.einsum','sign*eps/ws*np.einsum').replace('zeta/width*np.einsum','zeta/wd*np.einsum');exec(compile(src,str(f),'exec'),ns)
rows=[]
for grid,dt in [(1024,.0025),(2048,.00125)]:
 a=ns['run'](320,grid,dt,1);b=ns['run'](320,grid,dt,-1);R=-(a[0]-b[0])/(4*ns['eps']*ns['zeta']);q=np.arange(grid)*4/grid-2
 rows.append(dict(grid=grid,dt=dt,response=R,error_vs_ideal=R-ns['ideal'],max_energy_drift=max(abs(a[1]),abs(b[1])),source_area_before_normalization=float(np.sum(abs(q)<.075)*4/grid/.15),detector_area_before_normalization=float(np.sum(abs(q-.4)<.075)*4/grid/.15)))
print(json.dumps(rows,indent=2));p.with_name('results.json').write_text(json.dumps(rows,indent=2)+'\n');print('CP09_UNIT_AREA_CONTROL_EXECUTION_COMPLETE')

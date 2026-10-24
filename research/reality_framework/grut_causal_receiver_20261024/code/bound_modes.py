#!/usr/bin/env python3
"""Separate instantaneous bound modes as a diagnostic, not asymptotic particles."""
import json
import numpy as np
from scipy.linalg import eig_banded
from discrete_layer import ROOT,A0,A1,A2,A3


def project(name):
    run=json.loads((ROOT/f'{name}_RESULTS.json').read_text())
    table=np.load(ROOT/'data'/f'{name}_final.npz')['final'];r=table[:,0];h=run['parameters']['h'];n=len(r)
    u=r*(table[:,1]+1j*table[:,2]);p=r*(table[:,4]+1j*table[:,5])
    potential=table[:,3]**2+.1*(table[:,1]**2+table[:,2]**2)
    band=np.zeros((4,n));band[0]=potential-A0/h**2
    band[0,0]+=A2/h**2;band[0,-1]+=A2/h**2
    band[1,:-1]=-A1/h**2;band[1,0]+=A3/h**2;band[1,-2]+=A3/h**2
    band[2,:-2]=-A2/h**2;band[3,:-3]=-A3/h**2
    eigenvalues,vectors=eig_banded(band,lower=True,select='v',select_range=(0.,1.-1e-9),check_finite=True)
    U=vectors.T@u;P=vectors.T@p;weight=4*np.pi*h
    charge=-weight*np.imag(U.conj()*P);rows=[]
    # Independently apply the production seven-point stencil and odd ghosts,
    # rather than reusing the eigensolver's band storage for its residual.
    full=np.pad(vectors,((1,1),(0,0)))
    extended=np.concatenate([-full[3:0:-1],full,-full[-2:-5:-1]],axis=0)
    lap=A0*vectors
    for d,coefficient in [(1,A1),(2,A2),(3,A3)]:
        lap+=coefficient*(extended[4-d:n+4-d]+extended[4+d:n+4+d])
    applied=potential[:,None]*vectors-lap/h**2
    residual=float(np.max(abs(applied-vectors*eigenvalues[None,:]))) if len(eigenvalues) else 0.
    for j,value in enumerate(eigenvalues):
        omega=np.sqrt(value);core=float(sum(vectors[r<60,j]**2));q=float(charge[j])
        rows.append(dict(omega=float(omega),mass2_gap=float(1-value),canonical_Qbar=q,
            mode_norm_fraction_inside_60=core,mean_radius=float(np.dot(r,vectors[:,j]**2)),
            linear_mode_energy=float(.5*weight*(abs(P[j])**2+value*abs(U[j])**2))))
    allQ=float(-weight*np.sum(np.imag(u.conj()*p)))
    return dict(name=name,scope='Frozen instantaneous charged-field operator only. It includes the actual snapshot potential but is neither the full coupled fluctuation operator nor an asymptotic particle/bound-object classification.',
        bound_modes=rows,eigenvector_residual=residual,Qbar_physical_field=allQ,
        Qbar_in_instantaneous_subthreshold_modes=float(sum(charge)),
        Qbar_in_modes_with_99percent_norm_inside_60=float(sum(x['canonical_Qbar'] for x in rows if x['mode_norm_fraction_inside_60']>.99)),
        Qbar_in_other_subthreshold_modes=float(sum(x['canonical_Qbar'] for x in rows if x['mode_norm_fraction_inside_60']<=.99)))


def main():
    names=['P025_time','P0125','P025_remote_time','P025_remote']
    if (ROOT/'P00625_RESULTS.json').exists():names.append('P00625')
    rows=[project(name) for name in names]
    out=dict(scope='Late snapshot bound-component diagnostic. An instantaneous eigenvalue below m_X is not a future stable-object or particle certificate; time dependence, self-interaction and domain sensitivity remain.',rows=rows)
    (ROOT/'BOUND_MODE_DIAGNOSTICS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for row in rows:print(row['name'],'modes',len(row['bound_modes']),'Qbound',row['Qbar_in_instantaneous_subthreshold_modes'],
        'Qdeep',row['Qbar_in_modes_with_99percent_norm_inside_60'],'Qother',row['Qbar_in_other_subthreshold_modes'],flush=True)


if __name__=='__main__':main()

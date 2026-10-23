#!/usr/bin/env python3
"""Partial-wave rate table for one fixed, cold Q=200000 object."""
import json
from pathlib import Path
import time
import numpy as np
from scattering import Background, ROOT, MX, scatter, visible_spectral_probability

CM3S=.3893793721e-27*29979245800.


def smatrix(result):
    return np.array(result['S_real'])+1j*np.array(result['S_imag'])


def main():
    bg=Background();start=time.monotonic()
    momenta=np.unique(np.r_[np.geomspace(1e-5,.7071067811865476,20),
                             np.linspace(.1,.7071067811865476,13)])
    rows=[]
    for u in momenta:
        energy=np.sqrt(1+u*u)
        ellmax=max(3,int(np.ceil(u*6+6)))
        partials=[]
        for ell in range(ellmax+1):
            result,sol,waves=scatter(bg,energy,ell=ell,radius=30.,tol=1e-6,
                                     mesh=1001,return_solutions=True)
            prob=result['columns'][0]['probabilities'][0]
            sigmav=np.pi*(2*ell+1)*prob/(u*energy*MX*MX)*CM3S
            visible=visible_spectral_probability(bg,result,sol[0],waves,nr=2401,nq=301)
            photon=visible_spectral_probability(bg,result,sol[0],waves,nr=2401,nq=301,photon_coefficient=10.)
            partials.append(dict(ell=ell,P_conversion=prob,
                sigma_v_cm3_s=sigmav,P_visible_leptons=visible['probability'],
                P_visible_with_photon10=photon['probability'],
                flux_residual=result['columns'][0]['number_flux_residual'],
                BVP_residual=result['columns'][0]['max_BVP_residual']))
        total=sum(r['sigma_v_cm3_s'] for r in partials)
        row=dict(p_over_m=float(u),energy_over_m=float(energy),ellmax=ellmax,
            sigma_v_conversion_cm3_s=total,
            sigma_v_visible_leptons_cm3_s=float(sum(np.pi*(2*r['ell']+1)*r['P_visible_leptons']/(u*energy*MX*MX)*CM3S for r in partials)),
            sigma_v_visible_photon10_cm3_s=float(sum(np.pi*(2*r['ell']+1)*r['P_visible_with_photon10']/(u*energy*MX*MX)*CM3S for r in partials)),
            last_two_partial_fraction=sum(r['sigma_v_cm3_s'] for r in partials[-2:])/total,
            partials=partials)
        rows.append(row)
        out=dict(scope='Fixed cold Q=200000 background, linear waves, closed-s source band only; no hot-core lifetime or actual formation source.',
            Omega=bg.omega,profile_sha256=bg.hash,rows=rows,
            elapsed_seconds=time.monotonic()-start)
        (ROOT/'RATE_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print('rate',len(rows),float(u),total,flush=True)
    # Different radii, tolerances, seed meshes, profile precision, all open
    # channel columns and exact zero-potential controls.
    controls=[]
    for e in [1.001,np.sqrt(1.5)]:
        for radius,tol,mesh in [(25.,1e-6,801),(40.,2e-8,1201)]:
            for ell in [0,1,3]:
                r=scatter(bg,e,ell=ell,radius=radius,tol=tol,mesh=mesh,incoming=(0,1))
                S=smatrix(r)[:2]
                r['unitarity_max']=float(np.max(abs(S.conj().T@S-np.eye(2))))
                r['reciprocity_max']=float(np.max(abs(S-S.T)))
                controls.append(r)
    for ell in [0,1,4]:
        r=scatter(bg,np.sqrt(1.5),ell=ell,scale=0.,radius=30.,incoming=(0,1))
        S=smatrix(r)[:2]
        r['free_identity_error']=float(np.max(abs(S-np.eye(2))))
        controls.append(r)
    r=scatter(bg,1.8,radius=30.,tol=1e-7,incoming=(0,1,2))
    S=smatrix(r);r['unitarity_max']=float(np.max(abs(S.conj().T@S-np.eye(3))))
    r['reciprocity_max']=float(np.max(abs(S-S.T)));controls.append(r)
    r,sol,w=scatter(bg,np.sqrt(1.5),ell=0,radius=40.,tol=2e-8,mesh=1201,return_solutions=True)
    spectral=[visible_spectral_probability(bg,r,sol[0],w,nr=nr,nq=nq)
              for nr,nq in [(2401,301),(4801,601),(9601,1201)]]
    (ROOT/'SCATTERING_CONTROLS.json').write_text(json.dumps(dict(controls=controls,
           visible_kernel_quadrature=spectral,elapsed_seconds=time.monotonic()-start),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()

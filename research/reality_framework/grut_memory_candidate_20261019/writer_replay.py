#!/usr/bin/env python3
"""Independent exact-core / clock split and doubled-grid replay of M1-W."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

import autonomous_writer as original
from calculate import BASE, hamiltonian, state

HERE=Path(__file__).resolve().parent


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def rotate_strided(psi,bit,angle):
    stride=1<<bit
    view=psi.reshape(len(psi),-1,2,stride)
    a=view[:,:,0,:].copy();b=view[:,:,1,:].copy()
    theta=np.asarray(angle)
    if theta.ndim:
        theta=theta.reshape(len(psi),1,1)
    c,s=np.cos(theta),-1j*np.sin(theta)
    view[:,:,0,:]=c*a+s*b
    view[:,:,1,:]=c*b+s*a


def replay_scattering(points,dt,H,z,E,V):
    q=np.linspace(-16,16,points,endpoint=False);dq=q[1]-q[0]
    k=2*np.pi*np.fft.fftfreq(points,d=dq)
    kinetic=k*k/(2*original.M)
    packet=np.exp(-(q-original.Q0)**2/(4*original.SIGMA**2)+1j*original.P0*(q-original.Q0))
    packet/=np.linalg.norm(packet)
    psi=packet[:,None]*state(original.L,'readout',1)[None,:]
    gates=np.array([np.where(abs(q-i)<original.WIDTH,np.cos(np.pi*(q-i)/(2*original.WIDTH))**2,0.) for i in range(original.L)])
    ids=np.arange(len(H))
    def energy(w):
        fw=np.fft.fft(w,axis=0,norm='ortho')
        ek=float(np.sum(abs(fw)**2*kinetic[:,None]))
        es=float(np.sum(w.conj()*(w@H.T)).real)
        ev=0.
        for i in range(original.L):
            bit=original.L+i
            lo=ids[(ids&(1<<bit))==0];hi=lo|(1<<bit)
            ev+=2*original.OMEGA*float(np.sum(gates[i,:,None]*w[:,lo].conj()*w[:,hi]).real)
        return dict(clock=ek,spin=es,interaction=ev,total=ek+es+ev)
    initial=energy(psi)
    steps=int(round(original.TEND/dt));dt=original.TEND/steps
    Uhalf=(V*np.exp(-.5j*E*dt)[None,:])@V.T
    Ufull=(V*np.exp(-1j*E*dt)[None,:])@V.T
    kh=np.exp(-.5j*kinetic*dt);kf=kh*kh
    # Kinetic and core commute: combine their adjacent half-steps exactly.
    psi=np.fft.ifft(np.fft.fft(psi@Uhalf.T,axis=0,norm='ortho')*kh[:,None],axis=0,norm='ortho')
    for j in range(steps):
        for i in range(original.L):
            rotate_strided(psi,original.L+i,original.OMEGA*gates[i]*dt)
        u=Uhalf if j==steps-1 else Ufull
        kk=kh if j==steps-1 else kf
        psi=np.fft.ifft(np.fft.fft(psi@u.T,axis=0,norm='ortho')*kk[:,None],axis=0,norm='ortho')
        if (j+1)%100==0:
            print(json.dumps(dict(points=points,steps_complete=j+1,total_steps=steps)),flush=True)
    final=energy(psi)
    pq=np.sum(abs(psi)**2,axis=1)
    pk=np.sum(abs(np.fft.fft(psi,axis=0,norm='ortho'))**2,axis=1)
    rho=psi.T@psi.conj()
    return rho,dict(source_bit=-1,dt=dt,clock_points=points,initial_energy=initial,final_energy=final,
        total_energy_drift=final['total']-initial['total'],norm_error=abs(float(pq.sum())-1),
        spin_purity=float(np.trace(rho@rho).real),
        clock_gate_region_probability=float(pq[(q>-original.WIDTH-1)&(q<original.L-1+original.WIDTH+1)].sum()),
        clock_left_of_gates_probability=float(pq[q < -original.WIDTH].sum()),
        clock_edge_probability=float(pq[abs(q)>15].sum()),
        clock_high_momentum_probability=float(pk[abs(k)>.85*np.pi/dq].sum()),
        clock_mean_q=float(pq@q),clock_mean_p=float(pk@k),
        clock_negative_momentum_probability=float(pk[k<0].sum()))


def summarize(rho,minus,plusrho,plus,E,V,z,points,seconds):
    minus['storage']=original.storage_observables(rho,E,V,z)
    mem=.5+np.abs(np.array(plus['storage']['memory_means'])-np.array(minus['storage']['memory_means']))/4
    tw=np.linspace(0,2,101)
    ids=np.arange(len(E));Y=np.zeros((len(E),len(E)),dtype=complex)
    Y[ids,ids^2]=-1j*z[1]
    ye=V.T@Y@V;ze=V.T@(z[original.L,:,None]*V)
    rp=V.T@plusrho@V;rm=V.T@rho@V
    dE=E[:,None]-E[None,:]
    ys=[[float(np.sum(rr*ye.T*np.exp(-1j*dE*t)).real) for t in tw] for rr in [rp,rm]]
    read=.5+abs(np.array(ys[0])-ys[1])/4
    j=int(np.argmax(read));t=float(tw[j])
    pm=.5+abs(np.sum((rp-rm)*ze.T*np.exp(-1j*dE*t)).real)/4
    return dict(time_step=.0001,clock_points=points,plus=plus,minus=minus,
        memory_success=mem.tolist(),initial_all_cell_success=mem[:,0].tolist(),
        minimum_all_cell_success_through_t10=float(mem.min()),
        readout=dict(times=tw.tolist(),means=ys,success=read.tolist(),best_time=t,
                     best_success=float(read[j]),memory_at_best=float(pm)),wall_seconds=seconds)


def main():
    target=HERE/'writer_results_refined.json'
    if target.exists():raise SystemExit('Refusing to overwrite refined writer record')
    start=time.monotonic()
    old=json.loads((HERE/'writer_results_initial.json').read_text())
    coarse=old['runs'][0]
    H,z=hamiltonian(original.L,BASE);E,V=eigh(H)
    plusrho,plus=original.scattering(.0001,2048,1,H,z,E,V)
    plus['storage']=original.storage_observables(plusrho,E,V,z)
    result=dict(protocol_sha256=digest(HERE/'writer_protocol.json'),
        original_source_sha256=digest(HERE/'autonomous_writer.py'),
        source_sha256=digest(Path(__file__)),amendment_sha256=digest(HERE/'execution_amendment.json'),
        original_results_sha256=digest(HERE/'writer_results_initial.json'),
        runs=[coarse],scope=old['approximation'])
    for points in [2048,4096]:
        t0=time.monotonic()
        rho,minus=replay_scattering(points,.0001,H,z,E,V)
        run=summarize(rho,minus,plusrho,plus,E,V,z,points,time.monotonic()-t0)
        result['runs'].append(run)
        target.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(points=points,min_memory=run['minimum_all_cell_success_through_t10'],
            readout=run['readout']['best_success'],energy_drift=minus['total_energy_drift'],seconds=run['wall_seconds'])),flush=True)
    a,b,c=result['runs']
    def compare(x,y):
        return dict(max_memory_probability_difference=float(np.max(abs(np.array(x['memory_success'])-y['memory_success']))),
            max_readout_probability_difference=float(np.max(abs(np.array(x['readout']['success'])-y['readout']['success']))),
            energy_drift_difference=abs(x['minus']['total_energy_drift']-y['minus']['total_energy_drift']))
    result['refinement']=dict(time_and_algorithm=compare(a,b),spatial=compare(b,c),overall=compare(a,c))
    result['checks']=dict(write_and_hold=c['minimum_all_cell_success_through_t10']>=.9,
        all_cells_initial=min(c['initial_all_cell_success'])>=.9,
        readout=c['readout']['best_success']>.6 and c['readout']['memory_at_best']>=.9,
        norms=all(max(x['plus']['norm_error'],x['minus']['norm_error'])<1e-9 for x in result['runs']),
        probability_refinement=all(max(v['max_memory_probability_difference'],v['max_readout_probability_difference'])<=.01 for v in result['refinement'].values()),
        energy_refinement=all(v['energy_drift_difference']<=.1 for v in result['refinement'].values()))
    result['wall_seconds']=time.monotonic()-start
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=result['checks'],seconds=result['wall_seconds'])),flush=True)


if __name__=='__main__':main()

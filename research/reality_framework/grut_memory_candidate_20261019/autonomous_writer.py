#!/usr/bin/env python3
"""Time-independent quantum-particle writer, followed by stated scattering reduction.

The FFT box approximates R only during the bounded scattering calculation.
No transmission postselection, time-dependent physical drive, or reset is used.
Later storage neglects residual clock/gate overlap and reports that limitation.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from calculate import BASE, hamiltonian, state

HERE=Path(__file__).resolve().parent
L=4
M=.5
P0=100.
Q0=-4.
SIGMA=.35
WIDTH=.3
OMEGA=np.pi*(P0/M)/(2*WIDTH)
TEND=.06


def rotate(psi, low, high, angle):
    c,s=np.cos(angle),-1j*np.sin(angle)
    a=psi[:,low].copy();b=psi[:,high].copy()
    psi[:,low]=c*a+s*b
    psi[:,high]=c*b+s*a


def scattering(dt, points, source_bit, H, z, E, V):
    q=np.linspace(-16,16,points,endpoint=False)
    dq=q[1]-q[0]
    k=2*np.pi*np.fft.fftfreq(points,d=dq)
    kinetic=k*k/(2*M)
    packet=np.exp(-(q-Q0)**2/(4*SIGMA**2)+1j*P0*(q-Q0))
    packet/=np.linalg.norm(packet)
    spin=state(L,'readout',1)  # Same blank register/matter preparation for both inputs.
    psi=packet[:,None]*spin[None,:]
    ids=np.arange(len(H))
    pairs=[(ids[(ids&(1<<bit))==0],ids[(ids&(1<<bit))==0] | (1<<bit)) for bit in range(2*L)]
    gate=np.array([np.where(abs(q-i)<WIDTH,np.cos(np.pi*(q-i)/(2*WIDTH))**2,0.) for i in range(L)])
    diag=np.diag(H).copy()
    def energy_parts(wave):
        fw=np.fft.fft(wave,axis=0,norm='ortho')
        ek=float(np.sum(abs(fw)**2*kinetic[:,None]))
        es=float(np.sum(wave.conj()*(wave@H.T)).real)
        ev=0.
        if source_bit==-1:
            for i in range(L):
                lo,hi=pairs[L+i]
                ev+=2*OMEGA*float(np.sum(gate[i,:,None]*wave[:,lo].conj()*wave[:,hi]).real)
        return dict(clock=ek,spin=es,interaction=ev,total=ek+es+ev)
    initial_energy=energy_parts(psi)
    if source_bit==1:
        packet=np.fft.ifft(np.fft.fft(packet,norm='ortho')*np.exp(-1j*kinetic*TEND),norm='ortho')
        spin=V @ (np.exp(-1j*E*TEND)*(V.T@spin))
        psi=packet[:,None]*spin[None,:]
    else:
        steps=int(round(TEND/dt));dt=TEND/steps
        kinetic_half=np.exp(-.5j*kinetic*dt)
        # Symmetric local split for a half step of H_M1.
        step=dt/2
        diagonal_half=np.exp(-.5j*diag*step)
        memory_angles=[]
        for i in range(L):
            lo,_=pairs[L+i]
            F=BASE['epsilon']+(1-BASE['epsilon'])*(1-z[L+(i-1)%L,lo]*z[L+(i+1)%L,lo])/2
            memory_angles.append(-BASE['Delta']*F*step/2)
        def store_half():
            psi[:]*=diagonal_half[None,:]
            for i in range(L):
                rotate(psi,*pairs[i],BASE['hx']*step/2)
            for i in range(L):
                rotate(psi,*pairs[L+i],memory_angles[i][None,:])
            for i in reversed(range(L)):
                rotate(psi,*pairs[L+i],memory_angles[i][None,:])
            for i in range(L):
                rotate(psi,*pairs[i],BASE['hx']*step/2)
            psi[:]*=diagonal_half[None,:]
        for iteration in range(steps):
            psi[:]=np.fft.ifft(np.fft.fft(psi,axis=0,norm='ortho')*kinetic_half[:,None],axis=0,norm='ortho')
            store_half()
            for i in range(L):
                rotate(psi,*pairs[L+i],(OMEGA*gate[i]*dt)[:,None])
            store_half()
            psi[:]=np.fft.ifft(np.fft.fft(psi,axis=0,norm='ortho')*kinetic_half[:,None],axis=0,norm='ortho')
            if (iteration+1)%100==0:
                print(json.dumps(dict(dt=dt,bit=source_bit,steps_complete=iteration+1,total_steps=steps)),flush=True)
    final_energy=energy_parts(psi)
    pq=np.sum(abs(psi)**2,axis=1)
    pk=np.sum(abs(np.fft.fft(psi,axis=0,norm='ortho'))**2,axis=1)
    # All clock positions, including reflection, enter the trace.
    rho=psi.T@psi.conj()
    purity=float(np.trace(rho@rho).real)
    return rho,dict(source_bit=source_bit,dt=dt,clock_points=points,
                    initial_energy=initial_energy,final_energy=final_energy,
                    total_energy_drift=final_energy['total']-initial_energy['total'],
                    norm_error=abs(float(pq.sum())-1),spin_purity=purity,
                    clock_gate_region_probability=float(pq[(q>-WIDTH-1)&(q<L-1+WIDTH+1)].sum()),
                    clock_left_of_gates_probability=float(pq[q < -WIDTH].sum()),
                    clock_edge_probability=float(pq[abs(q)>15].sum()),
                    clock_high_momentum_probability=float(pk[abs(k)>.85*np.pi/dq].sum()),
                    clock_mean_q=float(pq@q),clock_mean_p=float(pk@k),
                    clock_negative_momentum_probability=float(pk[k<0].sum()))


def storage_observables(rho,E,V,z):
    ts=np.linspace(0,10,101)
    rhoE=V.T@rho@V
    dE=E[:,None]-E[None,:]
    def expectation(O):
        oe=V.T@O@V
        weights=rhoE*oe.T
        return [float(np.sum(weights*np.exp(-1j*dE*t)).real) for t in ts]
    memories=[expectation(np.diag(z[L+i])) for i in range(L)]
    ids=np.arange(len(E));Y=np.zeros((len(E),len(E)),dtype=complex)
    Y[ids,ids^2]=-1j*z[1]
    return dict(times=ts.tolist(),memory_means=memories,matter_Y=expectation(Y))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='writer_results_initial.json')
    parser.add_argument('--dt',nargs='+',type=float,default=[.0002,.0001])
    parser.add_argument('--points',type=int,default=2048)
    args=parser.parse_args()
    target=HERE/args.output
    if target.exists():raise SystemExit('Refusing to overwrite writer record')
    H,z=hamiltonian(L,BASE);E,V=eigh(H)
    start=time.monotonic()
    result=dict(protocol_sha256=hashlib.sha256((HERE/'writer_protocol.json').read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                parameters=dict(L=L,M=M,p0=P0,q0=Q0,sigma=SIGMA,w=WIDTH,Omega=OMEGA,T=TEND),
                runs=[],approximation='Finite FFT scattering calculation; subsequent spin storage neglects the reported small residual gate interaction, without transmission postselection. Not an infinite-time scattering certificate.')
    plusrho,plus=scattering(args.dt[0],args.points,1,H,z,E,V)
    plus['storage']=storage_observables(plusrho,E,V,z)
    for dt in args.dt:
        t0=time.monotonic()
        rho,minus=scattering(dt,args.points,-1,H,z,E,V)
        minus['storage']=storage_observables(rho,E,V,z)
        # Recompute readout on the registered .02 grid, rather than rounding storage .1 samples.
        tw=np.linspace(0,2,101)
        ids=np.arange(len(E));Y=np.zeros((len(E),len(E)),dtype=complex)
        Y[ids,ids^2]=-1j*z[1]
        ye=V.T@Y@V
        replus=V.T@plusrho@V;reminus=V.T@rho@V
        dE=E[:,None]-E[None,:]
        ys=[]
        for rr in [replus,reminus]:
            ys.append([float(np.sum(rr*ye.T*np.exp(-1j*dE*t)).real) for t in tw])
        read=.5+np.abs(np.array(ys[0])-ys[1])/4
        mem=.5+np.abs(np.array(plus['storage']['memory_means'])-np.array(minus['storage']['memory_means']))/4
        best=int(np.argmax(read));bestt=float(tw[best])
        ze=V.T@(z[L,:,None]*V)
        mem_at_best=.5+abs(np.sum((replus-reminus)*ze.T*np.exp(-1j*dE*bestt)).real)/4
        summary=dict(time_step=dt,clock_points=args.points,plus=plus,minus=minus,
                     memory_success=mem.tolist(),initial_all_cell_success=mem[:,0].tolist(),
                     minimum_all_cell_success_through_t10=float(mem.min()),
                     readout=dict(times=tw.tolist(),means=ys,success=read.tolist(),
                                  best_time=bestt,best_success=float(read[best]),memory_at_best=float(mem_at_best)),
                     wall_seconds=time.monotonic()-t0)
        result['runs'].append(summary)
        target.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(dt=dt,min_memory=float(mem.min()),readout=float(read[best]),
                 memory_at_read=float(mem_at_best),energy_drift=minus['total_energy_drift'],
                 residual_gate_probability=minus['clock_gate_region_probability'],seconds=summary['wall_seconds'])),flush=True)
    a,b=result['runs'][0],result['runs'][-1]
    result['refinement']=dict(max_memory_probability_difference=float(np.max(abs(np.array(a['memory_success'])-np.array(b['memory_success'])))),
        max_readout_probability_difference=float(np.max(abs(np.array(a['readout']['success'])-np.array(b['readout']['success'])))),
        energy_drift_difference=abs(a['minus']['total_energy_drift']-b['minus']['total_energy_drift']))
    result['checks']=dict(write_and_hold=b['minimum_all_cell_success_through_t10']>=.9,
        readout=b['readout']['best_success']>.6 and b['readout']['memory_at_best']>=.9,
        norms=all(max(x['plus']['norm_error'],x['minus']['norm_error'])<1e-9 for x in result['runs']),
        probability_refinement=max(result['refinement']['max_memory_probability_difference'],result['refinement']['max_readout_probability_difference'])<=.01,
        energy_refinement=result['refinement']['energy_drift_difference']<=.1)
    result['wall_seconds']=time.monotonic()-start
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=result['checks'],seconds=result['wall_seconds'])),flush=True)


if __name__=='__main__':main()

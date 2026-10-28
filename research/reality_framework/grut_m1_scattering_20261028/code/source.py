"""Measured outgoing monopole spectra from the prior fixed-M1 pulse traces."""
from common import *

def extract(tag,h,pad=4):
    path=WORK/'source'/tag/'state_trace.csv';histpath=WORK/'source'/tag/'state_history.csv'
    d=np.genfromtxt(path,names=True,delimiter=',');hist=np.genfromtxt(histpath,names=True,delimiter=',')
    dt=float(d['tau'][1]-d['tau'][0]);d=d[:-1];N=len(d);nf=pad*N;T=nf*dt
    z=d['x']+1j*d['y'];dz=d['dx']+1j*d['dy'];p=d['px']+1j*d['py']
    Z=dt*np.fft.fft(z,n=nf);DZ=dt*np.fft.fft(dz,n=nf);P=dt*np.fft.fft(p,n=nf)
    omega=2*np.pi*np.fft.fftfreq(nf,dt)
    # Actual C++ probe is the nearest cell center using positive std::lround.
    ip=int(np.floor(60/h));r=(ip+.5)*h;area=4*np.pi*r*r/G2
    charge=area/T*np.imag(np.conj(Z)*DZ)
    energy=-area/T*np.real(np.conj(P)*DZ)
    PS=dt*np.fft.fft(d['ps'],n=nf);DS=dt*np.fft.fft(d['ds'],n=nf)
    neutral_energy=-area/T*np.real(np.conj(PS)*DS)
    neutral_open=abs(omega)>np.sqrt(6)
    neutral_number=float(np.sum(neutral_energy[neutral_open]/abs(omega[neutral_open])))
    positive=omega < -1.;negative=omega > 1.;closed=abs(omega)<=1
    Np=float(charge[positive].sum());Nm=float(-charge[negative].sum())
    time_integral=float(area*np.sum(np.imag(np.conj(z)*dz))*dt)
    charge_fft=float(charge.sum())
    k=np.sqrt(omega[negative]**2-1);weight=-charge[negative]
    order=np.argsort(k);k=k[order];weight=weight[order]
    clipped=np.maximum(weight,0);cs=np.cumsum(clipped)
    quantiles=[float(np.interp(q*cs[-1],cs,k)) for q in (.1,.5,.9)]
    return {'tag':tag,'trace_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'history_sha256':hashlib.sha256(histpath.read_bytes()).hexdigest(),
            'sample_dt':dt,'padding':pad,'probe_radius':r,'Nplus_out':Np,'Nminus_out':Nm,
            'Nminus_over_initial_Q':Nm/QREF,'net_charge_out':Np-Nm,
            'Nphi_out_net_vacuum_reference':neutral_number,
            'phi_future_pair_number_after_one_vacuum_decay':neutral_number,
            'Nminus_plus_Nphi_over_Q_before_any_recapture':(Nm+neutral_number)/QREF,
            'phi_energy_out':float(neutral_energy[neutral_open].sum()),
            'subthreshold_signed_charge':float(charge[closed].sum()),
            'FFT_time_charge_identity_absolute':abs(charge_fft-time_integral),
            'net_out_vs_recorded_probe_relative':float(abs(charge_fft-hist['Qflux_probe'][-1]/G2)/max(1,abs(hist['Qflux_probe'][-1]/G2))),
            'recorded_receiver_charge':float(hist['Q_abs'][-1]/G2),
            'negative_weight_absolute_fraction':float(np.maximum(-weight,0).sum()/clipped.sum()),
            'k_quantiles_10_50_90':quantiles,'mean_negative_particle_energy':float(energy[negative].sum()/Nm),
            'above_neutral_opening_fraction':float(clipped[k>np.sqrt((np.sqrt(6)-OMEGA)**2-1)].sum()/clipped.sum())},k,weight

def main():
    out=[];arrays={}
    for tag,h in [('radial_main',.2),('radial_fine',.1),('pulse10',.2),('pulse20',.2),('pulse40_fine',.1)]:
        r,k,w=extract(tag,h);out.append(r);arrays[tag+'_k']=k;arrays[tag+'_negative_number_weight']=w
        print(json.dumps(r),flush=True)
    checks=[]
    for pad in (1,2,8):
        row,_,_=extract('radial_fine',.1,pad);checks.append(row)
    save('SOURCE_RESULTS.json',{'scope':__doc__,'rows':out,'padding_controls':checks,
                              'cosmological_preparation':False})
    np.savez_compressed(ROOT/'SOURCE_SPECTRA.npz',**arrays)

if __name__=='__main__':main()

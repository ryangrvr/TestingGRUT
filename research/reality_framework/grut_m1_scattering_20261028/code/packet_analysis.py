"""Directional massive-wave extraction and nonlinear/linear response comparison."""
from common import *
from packets import CASES
from scattering import scatter
from scipy.interpolate import PchipInterpolator

def directional(trace,h,dt,mass,pad=4,real_field=False):
    t=trace['tau'];spacing=float(t[1]-t[0]);N=len(t)-1;nf=pad*N;duration=nf*spacing
    r=(int(np.floor(40/h))+.5)*h
    if real_field:z=trace['s_minus1'][:-1];dz=trace['ds'][:-1]
    else:z=trace['x'][:-1]+1j*trace['y'][:-1];dz=trace['dx'][:-1]+1j*trace['dy'][:-1]
    Z=spacing*np.fft.fft(z,nf);DZ=spacing*np.fft.fft(dz,nf);w=2*np.pi*np.fft.fftfreq(nf,spacing)
    energy=2/dt*np.sin(abs(w)*dt/2);k2=energy**2-mass**2
    opened=(k2>0)&(h*h*k2<4)
    kj=np.zeros_like(w);kj[opened]=np.sqrt(k2[opened])*np.sqrt(1-h*h*k2[opened]/4)
    cosine=1-h*h*k2/2;U=r*Z
    derivative=(r*r-h*h)/r*DZ+cosine*Z
    plus=np.zeros_like(Z);minus=np.zeros_like(Z)
    plus[opened]=.5*(U[opened]+derivative[opened]/(1j*kj[opened]))
    minus[opened]=.5*(U[opened]-derivative[opened]/(1j*kj[opened]))
    factor=4*np.pi/G2/duration*kj
    if real_field:
        outgoing=np.where(w>0,abs(minus)**2,abs(plus)**2)*factor
        incoming=np.where(w>0,abs(plus)**2,abs(minus)**2)*factor
        return {'out':float(outgoing.sum()),'in':float(incoming.sum()),
                'energy_out':float(np.dot(outgoing,energy))}
    anti=w>0;particle=w<0
    values={'anti_in':float(np.sum(factor[anti]*abs(plus[anti])**2)),
            'anti_out':float(np.sum(factor[anti]*abs(minus[anti])**2)),
            'X_out':float(np.sum(factor[particle]*abs(plus[particle])**2)),
            'X_in':float(np.sum(factor[particle]*abs(minus[particle])**2))}
    ok=anti&opened
    return values,np.sqrt(k2[ok]),factor[ok]*abs(plus[ok])**2

def main():
    records={};arrays={}
    for tag in CASES:
        p=WORK/tag/'result.json'
        if not p.exists():raise RuntimeError('Missing packet '+tag)
        records[tag]=json.loads(p.read_text())
        arrays[tag+'_history']=np.genfromtxt(WORK/tag/'state_history.csv',names=True,delimiter=',')
        arrays[tag+'_trace']=np.genfromtxt(WORK/tag/'state_trace.csv',names=True,delimiter=',')
    out=[];response_tables={}
    for tag,record in records.items():
        if tag.startswith('baseline'):continue
        par=record['parameters'];h=par['h'];dt=par['dt'];trace=arrays[tag+'_trace'].copy()
        baseline='baseline005' if h==.05 else 'baseline01'
        if par['core']:
            b=arrays[baseline+'_trace']
            for n in trace.dtype.names:
                if n!='tau':trace[n]-=b[n]
        charged,k,weight=directional(trace,h,dt,1.)
        scalar=directional(trace,h,dt,np.sqrt(6),real_field=True)
        n=charged['anti_in'];coverage=(k>=.1)&(k<=3.3)
        # Build enough monopole samples to average the actual incoming packet.
        grid=np.unique(np.r_[np.linspace(.1,3.3,801),1.20743547399])
        if not par['core']:
            prediction=np.array([0.,float(weight[coverage].sum()/n),0.])
        else:
            if h not in response_tables:
                response=[]
                for kk in grid:
                    row=scatter(float(kk),h=h)
                    response.append([row[name] for name in ('P_X','P_antiparticle','P_phi')])
                response_tables[h]=np.array(response)
            response=response_tables[h]
            prediction=np.array([np.dot(weight[coverage],PchipInterpolator(grid,response[:,i],extrapolate=False)(k[coverage]))/n for i in range(3)])
        measured=np.array([charged['X_out'],charged['anti_out'],scalar['out']])/n
        hist=arrays[tag+'_history'];late=hist['tau']>=800
        if par['core']:
            bh=arrays[baseline+'_history'];delta=(hist['Q30']-bh['Q30'])/G2
            shift=float(delta[late].mean());spread=float(np.ptp(delta[late]))
        else:shift=None;spread=None
        row={'tag':tag,'h':h,'incoming_number_measured':n,'target_input_number':par['Nin'],
             'charged_channels':charged,'neutral_channels':scalar,
             'measured_P_X_anti_phi':measured.tolist(),
             'stationary_packet_weighted_P_X_anti_phi':prediction.tolist(),
             'absolute_response_differences':(measured-prediction).tolist(),
             'measured_number_flux_sum':float(measured.sum()),
             'late_core_charge_shift':shift,'late_core_charge_shift_range':spread,
             'charge_shift_from_outgoing_channels':float(-2*charged['X_out']-scalar['out']),
             'incident_number_covered_by_k_grid':float(weight[coverage].sum()/n)}
        out.append(row);print(json.dumps(row),flush=True)
    np.savez_compressed(ROOT/'PACKET_DATA.npz',**arrays)
    save('PACKET_RESULTS.json',{'scope':__doc__,'runs':records,'rows':out,
                              'actual_cosmological_capture_history':None})

if __name__=='__main__':main()

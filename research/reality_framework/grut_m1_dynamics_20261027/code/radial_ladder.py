"""Predeclared pulse-amplitude ladder, not a cosmological preparation scan."""
from radial import *

def main():
    cases={
      'pulse10':dict(h=.2,dt=.025,L=200.,ab=100.,amp=.1),
      'pulse20':dict(h=.2,dt=.025,L=200.,ab=100.,amp=.2),
      'pulse40':dict(h=.2,dt=.025,L=200.,ab=100.,amp=.4),
      'pulse40_fine':dict(h=.1,dt=.0125,L=200.,ab=100.,amp=.4),
    }
    out={};arrays={}
    for tag,params in cases.items():
        r=run_radial(tag,params)
        trace=np.genfromtxt(WORK/tag/'state_trace.csv',names=True,delimiter=',')
        r['max_probe_X_amplitude_over_v']=float(np.max(np.sqrt(trace['x']**2+trace['y']**2)))
        r['max_probe_phi_deviation_over_v']=float(np.max(abs(trace['s_minus1'])))
        r['lower_bound_emitted_negative_charge_per_Q']=max(0.,-r['absorbed_Q_over_original'])
        out[tag]=r;print(json.dumps(r),flush=True)
        arrays[tag+'_history']=np.genfromtxt(WORK/tag/'state_history.csv',names=True,delimiter=',')
        arrays[tag+'_final']=np.genfromtxt(WORK/tag/'state_final.csv',names=True,delimiter=',')
    save('EXCITATION_RESULTS.json',{'scope':__doc__,'cases':out,'physical_cosmic_p':None})
    np.savez_compressed(ROOT/'EXCITATION_DATA.npz',**arrays)

if __name__=='__main__':main()

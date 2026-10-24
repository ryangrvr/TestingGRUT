#!/usr/bin/env python3
"""Fresh saved-trace analysis; conditional linear mode numbers stay conditional."""
import hashlib,json
from pathlib import Path

import numpy as np
from scipy.signal.windows import tukey

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_core_conversion_20261023'
NAMES=['N025','N0125','N00625','N0125_far','N0125_box','N0125_abs']


def analyze(trace,parameters,alpha=.05,start=0.):
    tr=trace[trace[:,0]>=start];t=tr[:,0];dt=t[1]-t[0];n=len(t);T=n*dt
    r=round(parameters['probe']/parameters['h'])*parameters['h']
    window=tukey(n,alpha)
    value=tr[:,1]+1j*tr[:,2];derivative=tr[:,4]+1j*tr[:,5]
    momentum=tr[:,7]+1j*tr[:,8]
    F=dt*np.fft.fft(window*value);D=dt*np.fft.fft(window*derivative)
    P=dt*np.fft.fft(window*momentum);omega=-2*np.pi*np.fft.fftfreq(n,dt)
    area=4*np.pi*r*r
    energy=-area*np.real(P.conj()*D)/T;charge=area*np.imag(F.conj()*D)/T
    time_energy=-area*dt*np.sum(window**2*np.real(momentum.conj()*derivative))
    time_charge=area*dt*np.sum(window**2*np.imag(value.conj()*derivative))
    k2=omega**2-1;valid=k2>1e-4;k=np.sign(omega[valid])*np.sqrt(k2[valid])
    outgoing=.5*(F[valid]+(D[valid]+F[valid]/r)/(1j*k))
    incoming=F[valid]-outgoing
    eout=area*abs(omega[valid]*k)*abs(outgoing)**2/T
    ein=area*abs(omega[valid]*k)*abs(incoming)**2/T
    correction=abs((1+tr[:,3])**2-1+.1*abs(value)**2)
    controlled=k2[valid]>10*max(correction)
    number=eout/abs(omega[valid]);positive=omega[valid]>0
    hist_edges=np.r_[0.,np.geomspace(.01,4.,81),np.inf]
    hist=np.histogram(abs(k),hist_edges,weights=eout)[0]
    return dict(start=start,end=float(t[-1]),window_alpha=alpha,probe_radius=r,
        measured_windowed_net_energy=float(sum(energy)),measured_windowed_net_charge=float(sum(charge)),
        time_domain_windowed_energy=float(time_energy),time_domain_windowed_charge=float(time_charge),
        relative_energy_parseval_error=float(abs(sum(energy)-time_energy)/max(abs(time_energy),1e-100)),
        relative_charge_parseval_error=float(abs(sum(charge)-time_charge)/max(abs(time_charge),1e-100)),
        measured_mQ_over_E=float(sum(charge)/sum(energy)),
        max_local_mass2_correction=float(max(correction)),local_k2_cut=float(10*max(correction)),
        conditional_outgoing_linear_energy=float(sum(eout)),conditional_incoming_linear_energy=float(sum(ein)),
        conditional_incoming_over_outgoing=float(sum(ein)/sum(eout)),
        conditional_controlled_energy_fraction=float(sum(eout[controlled])/sum(eout)),
        conditional_all_band_antiparticle_number_fraction=float(sum(number[~positive])/sum(number)),
        conditional_controlled_antiparticle_number_fraction=float(sum(number[controlled&~positive])/max(sum(number[controlled]),1e-100)),
        endpoint_field_amplitude=float(abs(value[-1])),
        measured_propagating_net_energy=float(sum(energy[valid])),
        conditional_out_minus_in=float(sum(eout-ein)),
        spectrum_k_edges=[None if not np.isfinite(x) else float(x) for x in hist_edges],
        conditional_spectrum_outgoing_energy=[float(x) for x in hist])


def main():
    rows=[]
    for name in NAMES:
        run=json.loads((PARENT/f'{name}_RESULTS.json').read_text())
        path=PARENT/'data'/f'{name}_trace.npz';trace=np.load(path)['trace']
        sha=hashlib.sha256(trace.tobytes()).hexdigest();assert sha==run['trace_sha256']
        variants=[analyze(trace,run['parameters'],alpha) for alpha in [.025,.05,.1]]
        variants.append(analyze(trace,run['parameters'],.05,6000.))
        rows.append(dict(name=name,parameters=run['parameters'],trace_sha256=sha,
            trace_hash_matches_parent=True,trace_samples=len(trace),variants=variants,
            core60_late_mean_fraction=run['core60_late_mean_fraction']))
    first=np.loadtxt(PARENT/'data/N0125_box_final.csv',delimiter=',',skiprows=1)
    second=np.loadtxt(PARENT/'data/N0125_far_final.csv',delimiter=',',skiprows=1)
    common_fields_equal=np.array_equal(first,second)
    # Only the extraction radius changes; the physical history excludes probe flux columns.
    h1=np.loadtxt(PARENT/'data/N0125_box_history.csv',delimiter=',',skiprows=1)
    h2=np.loadtxt(PARENT/'data/N0125_far_history.csv',delimiter=',',skiprows=1)
    columns=[0,1,2,3,4,7,8,9,10,11,12,13,14]
    common_history_equal=np.array_equal(h1[:,columns],h2[:,columns])
    out=dict(scope='Reanalysis of six existing author production traces, not new formation simulations or independently certified asymptotic particles.',
        pair_caution='Signed current is measured; particle numbers require a free, outgoing, sufficiently long-duration interpretation. The all-band numbers do not meet this requirement uniformly.',
        rows=rows,matched_two_probe_control=dict(names=['N0125_box','N0125_far'],
            final_fields_bitwise_equal=bool(common_fields_equal),physical_histories_bitwise_equal=bool(common_history_equal),
            final_numeric_sha256=hashlib.sha256(first.tobytes()).hexdigest(),
            physical_history_numeric_sha256=hashlib.sha256(h1[:,columns].tobytes()).hexdigest()))
    (ROOT/'TRACE_AUDIT_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for row in rows:
        v=row['variants'][1]
        print(row['name'],'conditional controlled',v['conditional_controlled_energy_fraction'],
            'conditional anti',v['conditional_all_band_antiparticle_number_fraction'],flush=True)
    print('Two-probe control',out['matched_two_probe_control'],flush=True)


if __name__=='__main__':main()

"""Estimate full-field frequencies/growth, then compare with published predictions."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks

ROOT=Path(__file__).resolve().parent


def frequencies(t,signal,count,search_range):
    # No theoretical frequency is used to initialize this estimator.
    signal=signal/max(np.std(signal),1e-300)
    centered=signal-np.mean(signal)
    spectrum=abs(np.fft.rfft(centered*np.hanning(len(t))))**2
    bins=2*np.pi*np.fft.rfftfreq(len(t),d=t[1]-t[0])
    peaks,_=find_peaks(spectrum,distance=3)
    peaks=[i for i in peaks if search_range[0]<bins[i]<search_range[1]]
    peaks=sorted(peaks,key=lambda i:spectrum[i],reverse=True)[:count]
    if len(peaks)!=count:raise RuntimeError('Insufficient resolved frequency peaks')
    initial=np.sort(bins[peaks])
    centered_t=t-np.mean(t)
    def design(nu):
        return np.column_stack([np.ones(len(t)),*[f(w*centered_t) for w in nu for f in [np.cos,np.sin]]])
    def residual(nu):
        x=design(nu)
        coefficients=np.linalg.lstsq(x,signal,rcond=None)[0]
        return x@coefficients-signal
    spacing=2*np.pi/(t[-1]-t[0])
    fit=least_squares(residual,initial,bounds=(np.maximum(initial-spacing,search_range[0]),
                                             np.minimum(initial+spacing,search_range[1])),
                      ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=100)
    if not fit.success:raise RuntimeError('Frequency fit failed: '+fit.message)
    return dict(measured=np.sort(fit.x).tolist(),fft_initial=initial.tolist(),
                relative_RMS_residual=float(np.sqrt(np.mean(residual(fit.x)**2))))


def main():
    protocol=json.loads((ROOT/'protocol.json').read_text())
    manifest=json.loads((ROOT/'LOCAL_FREEZE.json').read_text())
    published=json.loads((ROOT.parent/'grut_theory_synthesis_20261020/DERIVATION_RESULTS.json').read_text())
    checks=[]
    def check(name,passed,detail):checks.append(dict(name=name,passed=bool(passed),detail=detail))
    for name,digest in manifest.items():
        check('Frozen '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,digest)
    output={};omega=.1**.25; f0=10**.25
    source_digest=hashlib.sha256((ROOT/'full_field.py').read_bytes()).hexdigest()
    for spec in protocol['runs']:
        data=json.loads((ROOT/'results'/(spec['id']+'.json')).read_text())
        trace=np.asarray(data['trace'])
        check(spec['id']+' numerical consistency',
              np.isfinite(trace).all() and data['source_sha256']==source_digest and
              data['maximum_relative_energy_drift']<protocol['analysis']['maximum_energy_drift'] and
              data['maximum_relative_charge_drift']<protocol['analysis']['maximum_charge_drift'],
              dict(energy_drift=data['maximum_relative_energy_drift'],charge_drift=data['maximum_relative_charge_drift']))
        row=dict(energy_drift=data['maximum_relative_energy_drift'],charge_drift=data['maximum_relative_charge_drift'])
        if spec['background']=='dense':
            k=spec['k']
            if spec['seed']=='charged':
                expected=[np.sqrt(k*k+3*omega*omega-sign*np.sqrt(9*omega**4+4*omega*omega*k*k)) for sign in [1,-1]]
                observed=frequencies(trace[:,0],trace[:,3],2,protocol['analysis']['frequency_search_range'])
            else:
                expected=[np.sqrt(k*k+f0*f0-1)]
                observed=frequencies(trace[:,0],trace[:,5],1,protocol['analysis']['frequency_search_range'])
            relative=abs(np.array(observed['measured'])/expected-1)
            row.update(expected=list(map(float,expected)),frequency_fit=observed,relative_errors=relative.tolist())
            check(spec['id']+' predicted frequency',max(relative)<protocol['analysis']['maximum_relative_frequency_error'],relative.tolist())
        elif spec['id'] in protocol['analysis']['growth_fit_windows']:
            start,end=protocol['analysis']['growth_fit_windows'][spec['id']]
            selected=(trace[:,0]>=start)&(trace[:,0]<=end)
            fit=np.polyfit(trace[selected,0],np.log(trace[selected,9]),1)
            # This scan is already published at 0cbb364, before these evolutions.
            expected=min(published['mixed_branch_prediction']['growth_scan'],key=lambda r:abs(r['k']-spec['k']))
            assert abs(expected['k']-spec['k'])<1e-12
            predicted=expected['growth']
            relative=abs(fit[0]/predicted-1)
            row.update(expected_growth=predicted,measured_growth=float(fit[0]),relative_error=float(relative),
                       fit_window=[start,end],largest_mode_norm=float(trace[:,9].max()),
                       largest_mode_norm_over_mean_F=float(np.max(trace[:,9]/trace[:,10])))
            check(spec['id']+' predicted growth',relative<protocol['analysis']['maximum_relative_growth_error'],float(relative))
        else:
            row.update(maximum_mode_norm_gain=float(trace[:,9].max()/trace[0,9]),
                       interpretation='Above the published instability band; report finite observed gain without inventing a post-run universal stability threshold.')
        output[spec['id']]=row
    comparisons={}
    base=output['dense_k02']['frequency_fit']['measured']
    for name in ['dense_k02_half_dt','dense_k02_double_grid']:
        comparisons[name]=dict(maximum_measured_frequency_difference=float(np.max(abs(np.array(base)-output[name]['frequency_fit']['measured']))))
    comparisons['growth_seed_amplitude']=dict(relative_change=float(abs(output['mixed_k04_large_seed']['measured_growth']/output['mixed_k04']['measured_growth']-1)))
    result=dict(experiment='Full nonlinear fields, plane-dependent periodic bulk, separate Fourier/RK4 implementation',
                completed=len(output),all_checks_pass=all(c['passed'] for c in checks),checks=checks,
                cases=output,comparisons=comparisons,
                scope='Author numerical tests of published model predictions. No finite-droplet, unrestricted-3D, quantum-lifetime or empirical verdict.')
    (ROOT/'SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()

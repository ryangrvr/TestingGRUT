"""Extend only the isolated quartic diagnostic; no full-system Floquet claim."""
from timing_common import *
from phase_and_orbit import orbit_diagnostic

def main():
    data=np.load(MODE/'MODE_PROFILES.npz');times=data['times']
    bg=data['b3_g0.5_bg'][np.argmin(abs(times-10))]
    ks=np.unique(np.r_[.01,.02,.05,.1,.2,.4,.6,2.13,np.linspace(.8,3.3,51),
                       np.arange(3.10,3.201,.005),np.arange(5.4,6.201,.01)])
    result=orbit_diagnostic(bg,ks)
    unstable=np.array([r['mu_conformal']>1e-6 for r in result['rows']]);indices=np.flatnonzero(unstable)
    groups=np.split(indices,np.flatnonzero(np.diff(indices)>1)+1);bands=[]
    for indices in groups:
        if not len(indices):continue
        lo=int(indices[0]);hi=int(indices[-1]);peak=max((result['rows'][i] for i in indices),key=lambda r:r['mu_conformal'])
        bands.append(dict(low_edge_bracket=[float(ks[max(0,lo-1)]),float(ks[lo])],
            high_edge_bracket=[float(ks[hi]),float(ks[min(len(ks)-1,hi+1)])],sampled_peak=peak,
            lower_edge_not_bracketed=bool(lo==0)))
    result['sampled_bands']=bands;result['positive_growth_threshold']=1e-6
    result['IR_mu_over_k']=[dict(k=r['k'],mu_over_k=r['mu_conformal']/r['k']) for r in result['rows'] if r['k']<=.05]
    result['scope']='Auxiliary isolated quartic orbit and supplied vacuum action only. New infrared/edge samples, not a Floquet theorem for the full reheating history or a selected cosmic seed.'
    save('QUARTIC_SCAN.json',result)
    assert len(bands)==2
    checks=dict(passed=result['symplectic_error']<1e-8 and result['circular_multiplier_error']<1e-8,
        symplectic=result['symplectic_error'],circular_multiplier_error=result['circular_multiplier_error'],
        sampled_bands=bands,IR_mu_over_k=result['IR_mu_over_k'])
    save('QUARTIC_CHECKS.json',checks);print((ROOT/'QUARTIC_CHECKS.json').read_text(),flush=True)

if __name__=='__main__':main()

"""Resolve the sampled low-k l4 resonance before a source-weighted rate claim."""
from common import *
from partials import run_grid
from scipy.interpolate import PchipInterpolator
from source import extract

def main():
    k=np.linspace(.30,.38,321)
    p,e=run_grid('resolved_l4_h005',k,lmax=8,h=.05)
    pf,ef=run_grid('resolved_l4_h0025',k,lmax=8,h=.025)
    factor=(2*np.arange(9)+1)[None,:]
    sigma=np.pi/k**2*np.sum(factor*p[:,:,0],axis=1)
    fine=np.pi/k**2*np.sum(factor*pf[:,:,0],axis=1)
    base=np.load(ROOT/'PARTIAL_WAVES.npz');coarsek=base['k'];coarse=base['sigma_X_over_mXminus2']
    dense=json.loads((ROOT/'PARTIAL_RESULTS.json').read_text())
    dk=np.array(dense['dense_grid_k']);ds=np.array(dense['dense_grid_sigma_X'])
    outside=((coarsek<k[0])|(coarsek>k[-1]))&((coarsek<dk[0])|(coarsek>dk[-1]))
    def interpolation(piece):
        kk=np.r_[coarsek[outside],k,dk];ss=np.r_[coarse[outside],piece,ds]
        order=np.argsort(kk);return PchipInterpolator(kk[order],ss[order])
    f=interpolation(sigma);ff=interpolation(fine);old=PchipInterpolator(coarsek,coarse)
    rows=[]
    for tag,h in [('radial_fine',.1),('pulse10',.2),('pulse20',.2),('pulse40_fine',.1)]:
        averages=[]
        for pad in (4,16,32):
            _,q,w=extract(tag,h,pad=pad);w=np.maximum(w,0);ok=(q>=coarsek[0])&(q<=coarsek[-1])
            v=q[ok]/np.sqrt(1+q[ok]**2)
            calc=lambda fun:float(np.dot(w[ok],fun(q[ok])*v)/w.sum())
            averages.append({'FFT_padding':pad,'sampled_grid_average':calc(old),'resolved_average':calc(f),
                              'resolved_resonance_finer_h_average':calc(ff),'covered_number_fraction':float(w[ok].sum()/w.sum())})
        rows.append({'source':tag,'controls':averages})
    np.savez_compressed(ROOT/'RESONANCE_DATA.npz',k=k,h005_probability=p,h0025_probability=pf,
                        h005_sigma=sigma,h0025_sigma=fine)
    out={'scope':__doc__,'k_interval':[.30,.38],'ell_range':[0,8],'rows':rows,
         'resonance_sigma_area_h005':float(np.trapezoid(sigma,x=k)),
         'resonance_sigma_area_h0025':float(np.trapezoid(fine,x=k)),
         'max_unitarity_error':float(max(e.max(),ef.max())),
         'all_possible_resonances_globally_certified':False}
    save('RESONANCE_CHECKS.json',out);print(json.dumps(out,indent=2))

if __name__=='__main__':main()

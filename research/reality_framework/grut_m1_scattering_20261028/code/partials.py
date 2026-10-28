"""Angular partial-wave conversion kernels and labelled source averages."""
from common import *
from scattering import scatter
from scipy.interpolate import PchipInterpolator

def run_grid(tag,ks,lmax=32,h=.05,Rfar=500.):
    path=WORK/(tag+'.npz')
    if path.exists():
        with np.load(path) as a:return a['probabilities'],a['unitarity']
    probs=np.zeros((len(ks),lmax+1,3));errors=np.zeros((len(ks),lmax+1))
    for i,k in enumerate(ks):
        for ell in range(lmax+1):
            r=scatter(float(k),ell=ell,h=h,Rfar=Rfar)
            probs[i,ell]=[r['P_X'],r['P_antiparticle'],r['P_phi']];errors[i,ell]=r['unitarity_error']
        if (i+1)%40==0:print(json.dumps({'grid':tag,'completed_k':i+1,'total':len(ks)}),flush=True)
    np.savez_compressed(path,k=ks,probabilities=probs,unitarity=errors)
    return probs,errors

def main():
    ks=np.unique(np.r_[np.linspace(.1,1.2,111),np.linspace(1.21,3.3,140),[.4,.7,1.,1.25,1.5,2.,3.]])
    probs,err=run_grid('partial_grid',ks)
    ell=np.arange(probs.shape[1]);weight=(2*ell+1)[None,:]
    sx=np.pi/ks**2*np.sum(weight*probs[:,:,0],axis=1)
    sphi=np.pi/ks**2*np.sum(weight*probs[:,:,2],axis=1)
    arrays={'k':ks,'P_X':probs[:,:,0],'P_antiparticle':probs[:,:,1],'P_phi':probs[:,:,2],
            'sigma_X_over_mXminus2':sx,'sigma_phi_over_mXminus2':sphi,'unitarity':err}
    np.savez_compressed(ROOT/'PARTIAL_WAVES.npz',**arrays)
    table=[]
    for k in (.2,.4,.7,1.,1.15,1.25,1.5,2.,3.):
        i=int(np.argmin(abs(ks-k)))
        tails={str(c):float(np.sum(weight[0,c+1:]*probs[i,c+1:,0])/max(np.sum(weight[0]*probs[i,:,0]),1e-300)) for c in (8,16,24)}
        table.append({'k':float(ks[i]),'sigma_charge_exchange_mXminus2':float(sx[i]),
                       'sigma_neutral_mXminus2':float(sphi[i]),'charge_exchange_tail_fractions':tails})
    source=np.load(ROOT/'SOURCE_SPECTRA.npz');averages=[]
    for name in ('radial_fine','pulse10','pulse20','pulse40_fine'):
        k=source[name+'_k'];w=source[name+'_negative_number_weight'];positive=np.maximum(w,0)
        ok=(k>=ks[0])&(k<=ks[-1]);velocity=k[ok]/np.sqrt(1+k[ok]**2)
        cross=PchipInterpolator(ks,sx)(k[ok]);neutral=PchipInterpolator(ks,sphi)(k[ok])
        averages.append({'source':name,'calculated_number_fraction':float(positive[ok].sum()/positive.sum()),
                         'sigma_v_X_covered_mXminus2':float(np.dot(positive[ok],cross*velocity)/positive.sum()),
                         'sigma_v_phi_covered_mXminus2':float(np.dot(positive[ok],neutral*velocity)/positive.sum()),
                         'single_collision_plane_wave_cross_section_not_a_global_optical_depth':True})
    # Denser source-dominated interval checks quadrature/interpolation and angular tail.
    dense=np.linspace(.65,1.15,101);denseprob,denseerr=run_grid('source_dense_grid',dense)
    dense_sx=np.pi/dense**2*np.sum(weight*denseprob[:,:,0],axis=1)
    checks=[]
    for k in (.4,1.,1.5,3.):
        p48,_=run_grid('tail48_'+str(k),np.array([k]),lmax=48)
        p25,_=run_grid('fine025_'+str(k),np.array([k]),lmax=32,h=.025)
        pfar,_=run_grid('far1000_'+str(k),np.array([k]),lmax=32,Rfar=1000.)
        calc=lambda p:float(np.pi/k**2*np.sum((2*np.arange(p.shape[1])+1)*p[0,:,0]))
        base=float(PchipInterpolator(ks,sx)(k))
        checks.append({'k':k,'base_sigma':base,'l48_sigma':calc(p48),'h025_sigma':calc(p25),'far1000_sigma':calc(pfar)})
    save('PARTIAL_RESULTS.json',{'scope':__doc__,'lmax':32,'h':.05,'momentum_grid_count':len(ks),
         'max_unitarity_error':float(err.max()),'table':table,'source_averages':averages,
         'refinements':checks,'dense_grid_k':dense.tolist(),'dense_grid_sigma_X':dense_sx.tolist(),
         'dense_grid_max_unitarity_error':float(denseerr.max()),'actual_cosmological_capture_history':None})
    print(json.dumps({'table':table,'source_averages':averages,'refinements':checks},indent=2))

if __name__=='__main__':main()

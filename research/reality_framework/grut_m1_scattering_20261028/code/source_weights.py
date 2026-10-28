"""Source-weighted angular kernels and a denser momentum-quadrature control."""
from common import *
from scipy.interpolate import PchipInterpolator

def main():
    result=json.loads((ROOT/'PARTIAL_RESULTS.json').read_text());src=np.load(ROOT/'SOURCE_SPECTRA.npz')
    with np.load(ROOT/'PARTIAL_WAVES.npz') as d:k=d['k'];sigma=d['sigma_X_over_mXminus2']
    densek=np.array(result['dense_grid_k']);densesigma=np.array(result['dense_grid_sigma_X'])
    mask=(k<densek[0])|(k>densek[-1]);kk=np.r_[k[mask],densek];ss=np.r_[sigma[mask],densesigma]
    order=np.argsort(kk);fine=PchipInterpolator(kk[order],ss[order]);base=PchipInterpolator(k,sigma)
    out=[]
    for tag in ('radial_fine','pulse10','pulse20','pulse40_fine'):
        q=src[tag+'_k'];w=np.maximum(src[tag+'_negative_number_weight'],0);ok=(q>=k[0])&(q<=k[-1]);v=q[ok]/np.sqrt(1+q[ok]**2)
        a=float(np.dot(w[ok],base(q[ok])*v)/w.sum());b=float(np.dot(w[ok],fine(q[ok])*v)/w.sum())
        out.append({'source':tag,'base_sigma_v_mXminus2_covered':a,'dense_sigma_v_mXminus2_covered':b,
                    'relative_quadrature_change':b/a-1,'covered_number_fraction':float(w[ok].sum()/w.sum())})
    save('SOURCE_WEIGHT_CHECKS.json',{'rows':out,'uncomputed_momentum_tails_not_extrapolated':True})
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()

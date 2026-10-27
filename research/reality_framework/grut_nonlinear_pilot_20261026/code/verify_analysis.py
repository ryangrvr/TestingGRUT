"""Tests of torus topology, energy Parseval and reference charge bookkeeping."""
from analyze import *
from negative_candidates import negative_catalog

def main():
    groups={};N=24;L=24.;mask=np.zeros((N,N,N),bool);mask[:2,5:8,5:8]=True;mask[-2:,5:8,5:8]=True
    lab,n=periodic_labels(mask)
    groups['periodic_component_merge']=dict(passed=n==1 and len(np.unique(lab[mask]))==1,count=n)
    rng=np.random.default_rng(20261027);q=rng.normal(size=(3,N,N,N));p=rng.normal(size=q.shape);p[1]+=2*q[0];q[2]+=C.NU
    sp=spectra(q,p,L);ee=local_energy(q,p,L).mean();ef=energy(q,p,L)
    groups['energy_and_spectral_bookkeeping']=dict(passed=abs(ee/ef-1)<1e-12 and sp['charge_Parseval_error']<1e-10 and sp['vacuum_reference']['charge_identity_relative']<1e-12,
        energy_relative=abs(ee/ef-1),charge_Parseval=sp['charge_Parseval_error'],reference_charge=sp['vacuum_reference']['charge_identity_relative'],
        scope='Reference particle charge identity holds even in interacting fields and does not certify particle numbers.')
    xx=np.fft.fftfreq(N)*N;rr=np.sqrt(xx[:,None,None]**2+xx[None,:,None]**2+xx[None,None,:]**2)
    q=np.zeros((3,N,N,N));q[0]=12*np.exp(-rr**2/18);q[2]=C.NU*(1-.95*np.exp(-rr**2/32));p=np.zeros_like(q);p[1]=.85*q[0]
    cat,lab=catalog(q,p,L,threshold=4.,smooth=.5)
    best=cat['candidates'][0];dist=np.linalg.norm(np.minimum(best['centre'],L-np.array(best['centre'])))
    groups['boundary_centred_candidate']=dict(passed=cat['candidate_count']==1 and dist<1e-8 and best['compact'],count=cat['candidate_count'],centre_error=dist,
        scope='Synthetic Gaussian candidate tests wrapping and moments, not a stationary soliton.')
    minus=p.copy();minus[1]*=-1;neg=negative_catalog(q,minus,L,threshold=4.)
    # The negative routine uses smoothing length1, so compare at that same
    # length rather than the earlier .5 boundary-moment check.
    positive,_=catalog(q,p,L,threshold=4.,smooth=1.)
    ce=abs(neg['candidates'][0]['opposite_Q_mask_magnitude']/positive['candidates'][0]['Q_mask']-1)
    groups['boundary_centred_candidate']['charge_conjugation_error']=ce
    groups['boundary_centred_candidate']['passed']=bool(groups['boundary_centred_candidate']['passed'] and ce<1e-12)
    save('ANALYSIS_CHECKS.json',dict(groups=groups))
    print(groups,flush=True)
    if not all(g['passed'] for g in groups.values()):raise SystemExit(1)

if __name__=='__main__':main()

"""Time-averaged field exports and a labelled classical frequency action."""
from pt_common import *
from pt_production import CASES,NAMES

def legacy_averages():
    z=np.load(PARENT/'FIELD_CURVES.npz');rows=[]
    tags=sorted(k.removesuffix('_tau') for k in z.files if k.endswith('_tau'))
    for tag in tags:
        t=z[tag+'_tau'];s=z[tag+'_eta'];sel=(t>=700-1e-9)&(t<=800+1e-9)
        row=dict(tag=tag,tau_window=[700.,800.],samples=int(sel.sum()),statistics={})
        for name in ['energy_outside_spatial_mean_fraction','phi_plus_interaction_energy_fraction','spatial_mean_charge_fraction']:
            v=z[tag+'_'+name];values=v[sel];tt=t[sel];ss=s[sel]
            row['statistics'][name]=dict(cosmic_time_mean=float(np.trapezoid(values,tt)/(tt[-1]-tt[0])),
                conformal_time_mean=float(np.trapezoid(values,ss)/(ss[-1]-ss[0])),minimum=float(values.min()),maximum=float(values.max()),last=float(values[-1]))
        rows.append(row)
    save('LEGACY_ORBIT_AVERAGES.json',dict(rows=rows,scope='Original seven public curve arrays, including their stated decimation. Envelopes are sampled, not certified continuous extrema. No old trajectory is rerun.'))

def peak(w,P,sign):
    mask=w*sign>0;ww=w[mask];pp=P[mask]
    if not np.any(pp>0):return None
    i=int(np.argmax(pp));lo=hi=i
    while lo>0 and pp[lo]>pp[i]/2:lo-=1
    while hi<len(pp)-1 and pp[hi]>pp[i]/2:hi+=1
    return dict(omega=float(ww[i]),FWHM=float(abs(ww[hi]-ww[lo])),mean_abs_omega=float(np.sum(abs(ww)*pp)/np.sum(pp)))

def window(tag,lo,hi):
    directory=WORK/tag;info=read(f'RUN_{tag}.json');meta=np.load(directory/'modes.npz');tr=np.load(directory/'trace.npz')
    t=meta['s'];sel=(t>=lo-1e-9)&(t<=hi+1e-9);ix=np.flatnonzero(sel);psi=np.load(directory/'psi.npy',mmap_mode='r')
    samples=len(ix);ds=t[1]-t[0];w=2*np.pi*np.fft.fftfreq(samples,ds);h=np.hanning(samples);norm=samples*np.sum(h*h)
    weight=h*h/np.sum(h*h);N,L=int(meta['N']),float(meta['L']);Q=float(meta['nQ']);kc=meta['kcont'];kl=meta['kfd2']
    sid=np.floor(kc/(2*np.pi/L)+.5).astype(int);shells=[];spectra=[];totals=np.zeros(4)
    for shell in np.unique(sid):
        cols=np.flatnonzero(sid==shell);P=np.zeros(samples)
        # Bound temporary memory even in the largest volume.
        for start in range(0,len(cols),512):
            a=np.asarray(psi[np.ix_(ix,cols[start:start+512])],dtype=np.complex128)
            X=fft(h[:,None]*a,axis=0);P+=np.sum(abs(X)**2,axis=1)/norm
        plus=float(np.sum((w*P)[w>0])/Q);minus=float(np.sum((-w*P)[w<0])/Q)
        row=dict(shell=int(shell),k_center=float(shell*2*np.pi/L),k_min=float(kc[cols].min()),k_max=float(kc[cols].max()),
            kfd2_mean=float(kl[cols].mean()),modes=len(cols),negative_frequency_action_per_net=minus,positive_frequency_action_per_net=plus,
            positive_peak=peak(w,P,1),negative_peak=peak(w,P,-1))
        shells.append(row);spectra.append(P)
        if shell==0:totals[0]+=minus;totals[1]+=plus
        else:totals[2]+=minus;totals[3]+=plus
    vals=tr['values'][ix];names=list(tr['names']);stats={}
    for j,name in enumerate(names):
        v=vals[:,j];stats[name]=dict(hann2_mean=float(weight@v),mean=float(np.mean(v)),minimum=float(v.min()),maximum=float(v.max()))
    zminus,zplus,wminus,wplus=totals
    expectedQ=float(weight@vals[:,names.index('nQ')])/Q
    result=dict(tag=tag,window_s=[float(t[ix[0]]),float(t[ix[-1]])],window_tau=[float(clock(t[ix[0]])),float(clock(t[ix[-1]]))],
        frequency_resolution=2*np.pi/(samples*ds),nQ=Q,coherent_negative_frequency_action=zminus,
        wave_negative_frequency_action=wminus,total_negative_frequency_action=zminus+wminus,
        spectral_charge_ratio=zplus+wplus-zminus-wminus,equal_time_full_charge_ratio=expectedQ,
        omitted_k_gt8_Hartree_negative_reference=stats['tail_H']['hann2_mean'],
        low_wave_fraction_kcenters_le1p6=sum(x['negative_frequency_action_per_net'] for x in shells if x['shell'] and x['k_center']<=1.6)/wminus,
        fundamental_shell_wave_fraction=sum(x['negative_frequency_action_per_net'] for x in shells if x['shell']==1)/wminus,
        equal_time=stats,shells=shells,
        interpretation='Hann-weighted |frequency| classical field action, coherent zero mode separate. No free-particle basis, quantum population or abundance certification.',physical_nbar=None)
    key=f'{tag}_{lo:g}_{hi:g}';save('WINDOW_'+key+'.json',result)
    np.savez_compressed(ROOT/('SPECTRUM_'+key+'.npz'),omega=w,power_by_shell=np.array(spectra),shell=np.array([x['shell'] for x in shells]),Q=Q)
    print('WINDOW',key,'coherent',zminus,'waves',wminus,'total',zminus+wminus,'charge',result['spectral_charge_ratio'],flush=True)
    return result

def main():
    legacy_averages();tags=sys.argv[1:] or [c['tag'] for c in CASES]
    rows=[]
    for tag in tags:
        wins=[(14.,24.),(24.,34.)]
        if tag=='base48':wins=[(0.,6.),(6.,12.),*wins]
        for lo,hi in wins:rows.append(window(tag,lo,hi))
    batch=dict(windows=[dict(tag=r['tag'],s=r['window_s'],negative_action=r['total_negative_frequency_action']) for r in rows])
    (WORK/('analysis_'+hashlib.sha256(' '.join(tags).encode()).hexdigest()[:12]+'.json')).write_text(json.dumps(batch,indent=2)+'\n')

if __name__=='__main__':main()

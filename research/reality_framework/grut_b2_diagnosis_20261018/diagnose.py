"""B2 unchanged-law diagnosis. All physical propagation retains all states.

The pinned B2 builder/propagator are reused and disclosed, not independent
implementations. Spectral pinching and isolated cells are diagnostic references.
"""
import argparse, hashlib, importlib.util, json, math, pathlib, resource, time
from fractions import Fraction as F
import numpy as np
from scipy.linalg import norm
from scipy.sparse import csr_matrix

HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=HERE.parent/'grut_b2_collision_channel_20261018'
spec=importlib.util.spec_from_file_location('pinned_b2',(OLD/'checks.py'))
b2=importlib.util.module_from_spec(spec);spec.loader.exec_module(b2)
b1=b2.b1

def pair(z):return [float(np.real(z)),float(np.imag(z))]
def arr(x):return [float(v) for v in x]


def physical_account():
    model=b1.PathModel(8,sites=(0,7),resource=2)
    source=int(np.argmax(abs(model.initial())))
    prop=b2.Propagator(model.H,source)
    theta=np.arange(193)/32;times=theta/b1.JPAIR
    states=prop.v@(np.exp(-1j*prop.e[:,None]*times)*prop.coeff[:,None])
    probs=abs(states)**2
    q=model.pairs;f=np.array([ff for p,m,ff in model.states])
    psets=[set(p) for p,m,ff in model.states];msets=[set(m) for p,m,ff in model.states]
    co=np.array([p==m for p,m,ff in model.states])
    adj=np.array([len(p)==2 and p==m and p[1]==p[0]+1 for p,m,ff in model.states])
    R=np.array([ff==0 and p==m and len(p)==2 and p[0] in (1,2) and p[1] in (5,6)
                for p,m,ff in model.states])
    P=(q==0)|((q==1)&co)|((q==2)&co&~adj)
    assert int(np.sum(P))==38
    Lm=np.array([(int(0 in p)+int(0 in m))/2 for p,m in zip(psets,msets)])
    Rm=np.array([(int(7 in p)+int(7 in m))/2 for p,m in zip(psets,msets)])
    FL=(f&1)>0;FR=(f&2)>0;LR=Lm+FL;RR=Rm+FR
    fields=np.array([sum(e*e for e in ef) for ef in model.fields])
    def prob(mask):return np.sum(probs[mask],axis=0)
    def mean(d):return d@probs
    pq=[prob(q==a) for a in range(3)]
    hstates=model.H@states
    masks=[(q==a)[:,None] for a in range(3)]
    z01=np.sum((states*masks[1]).conj()*(model.C@(states*masks[0])),axis=0)
    z12=np.sum((states*masks[2]).conj()*(model.C@(states*masks[1])),axis=0)
    i01=2*z01.imag/b1.JPAIR;i12=2*z12.imag/b1.JPAIR
    dp=[2*np.sum((states*m).conj()*hstates,axis=0).imag/b1.JPAIR for m in masks]
    assert np.max(abs(dp[0]+i01))<2e-9
    assert np.max(abs(dp[1]-i01+i12))<2e-9
    assert np.max(abs(dp[2]-i12))<2e-9
    sitecurr=[]
    cc=model.C.tocoo()
    for bit,number in [(1,FL),(2,FR)]:
        keep=(f[cc.row]^f[cc.col])==bit
        ce=csr_matrix((cc.data[keep],(cc.row[keep],cc.col[keep])),shape=model.H.shape)
        off=(~number)[:,None];on=number[:,None]
        z=np.sum((states*off).conj()*(ce@(states*on)),axis=0)
        sitecurr.append(2*z.imag/b1.JPAIR)
    hout=model.T@states
    joutL=-2*np.sum((LR[:,None]*states).conj()*hout,axis=0).imag/b1.JPAIR
    joutR=-2*np.sum((RR[:,None]*states).conj()*hout,axis=0).imag/b1.JPAIR
    # Conversion cannot change either endpoint's total resource.
    assert np.max(abs(((LR[:,None]-LR[None,:])*model.C.toarray())),initial=0)==0
    assert np.max(abs(((RR[:,None]-RR[None,:])*model.C.toarray())),initial=0)==0
    amp=states[source]*np.exp(1j*(2*b1.OMEGA)*times)
    columns=['theta','p_q0','p_q1','p_q2','p_fuel11','p_fuel10','p_fuel01','p_fuel00',
             'fuel_left','fuel_right','p_q1_colocated','p_q2_both_colocated','p_adjacent_copairs',
             'p_output_R','left_total_resource','right_total_resource','left_matter_resource',
             'right_matter_resource','I01_dtheta','I12_dtheta','z01_real','z01_imag',
             'z12_real','z12_imag','left_creation_dtheta','right_creation_dtheta',
             'left_resource_out_dtheta','right_resource_out_dtheta','return_amplitude_rel_real',
             'return_amplitude_rel_imag','bare_nonP_probability','mean_electric_energy']
    values=[theta,*pq,prob(f==3),prob(f==1),prob(f==2),prob(f==0),mean(FL),mean(FR),
            prob((q==1)&co),prob((q==2)&co),prob(adj),prob(R),mean(LR),mean(RR),mean(Lm),mean(Rm),
            i01,i12,z01.real,z01.imag,z12.real,z12.imag,*sitecurr,joutL,joutR,amp.real,amp.imag,
            prob(~P),mean(fields)]
    table=np.array(values).T
    assert table.shape==(193,len(columns))
    assert np.max(abs(np.sum(probs,axis=0)-1))<1e-10
    assert np.max(abs(np.sum(pq,axis=0)-1))<1e-10
    integer=[dict(zip(columns,arr(table[k]))) for k in range(0,193,32)]
    start=64;end=192
    balances=[]
    for label,fd,md,ld in [('left',mean(FL),mean(Lm),mean(LR)),('right',mean(FR),mean(Rm),mean(RR))]:
        df=fd[end]-fd[start];dm=md[end]-md[start];dl=ld[end]-ld[start]
        assert abs(df-dl+dm)<1e-12
        balances.append({'source':label,'fuel_gain_2_to_6':float(df),'matter_change_at_source':float(dm),
                         'net_resource_returned_across_bond':float(dl),
                         'initial_source_matter_resource_at_2':float(md[start]),
                         'final_source_matter_resource_at_6':float(md[end]),
                         'identity_defect':float(df-dl+dm)})
    perinterval=[]
    for a in range(6):
        aa=32*a;bb=aa+32
        perinterval.append({'theta_interval':[a,a+1],'net_creation_0_to_1':float(pq[0][aa]-pq[0][bb]),
            'net_creation_1_to_2':float(pq[2][bb]-pq[2][aa]),
            'left_resource_inflow':float(mean(LR)[bb]-mean(LR)[aa]),
            'right_resource_inflow':float(mean(RR)[bb]-mean(RR)[aa])})
    isolated=[]
    for t in range(7):
        pp=576/745*math.sin(math.sqrt(745)*t/24)**2
        isolated.append({'theta':t,'per_cell_matter':pp,'p_q0':(1-pp)**2,
                         'p_q1':2*pp*(1-pp),'p_q2':pp**2,'endpoint_resource_each':1})
    summary={'sector_definition':'q counts matter-pair quanta; co-location is separately reported',
             'integer_time_points':integer,'source_balance_2_to_6':balances,
             'signed_current_integrals_from_exact_continuity':perinterval,
             'isolated_cell_mathematical_reference':isolated,
             'state_error_account_at_endpoint':prop.account(),
             'failed_B2_endpoint_unchanged':True}
    return model,prop,states,table,columns,summary,{'q0':q==0,'q1':q==1,'q2':q==2,'R':R,'bareP':P}


def spectral_account(model,prop,states,masks):
    A=F(13,46080);J=F(1,3840);Om=2-2*A;tau2=F(5,48)
    high_edge=F(9,2)-A-tau2-2*J;gap=high_edge-2*Om
    wbound=2*J*J/(gap*gap)
    assert wbound==F(32,37027225)
    high=prop.e>=float(high_edge);assert np.sum(~high)==38
    weights=abs(prop.coeff)**2;w=float(np.sum(weights[high]));assert w<float(wbound)
    delta=prop.orth
    eta=prop.rho/math.sqrt(1-delta)+2*math.sqrt(1+delta)*max(abs(prop.e))*(1/math.sqrt(1-delta)-1)+prop.repr
    groups=[];group=[0]
    for i in range(1,len(prop.e)):
        if prop.e[i]-prop.e[i-1]<=2*eta:group.append(i)
        else:groups.append(group);group=[i]
    groups.append(group)
    T=6/b1.JPAIR;phase=np.exp(-1j*prop.e*T)
    pinched=np.zeros(model.dim);group_weights=[];low_groups=[];high_return=0j
    for ids in groups:
        vec=prop.v[:,ids]@(phase[ids]*prop.coeff[ids])
        pinched+=abs(vec)**2
        group_weights.append(float(np.sum(weights[ids])))
        contribution=sum(weights[ids]*np.exp(-1j*(prop.e[ids]-2*b1.OMEGA)*T))
        if not high[ids[0]]:
            assert not np.any(high[ids])
            low_groups.append({'energy_interval':[float(prop.e[ids[0]]-eta),float(prop.e[ids[-1]]+eta)],
                               'dimension':len(ids),'weight':group_weights[-1],
                               'return_amplitude_rel':pair(contribution)})
        else:high_return+=contribution
    assert abs(np.sum(pinched)-1)<1e-10
    probs=abs(states[:,-1])**2
    comparisons={name:{'physical_endpoint':float(np.sum(probs[mask])),
                       'spectral_block_pinched_reference':float(np.sum(pinched[mask])),
                       'resolved_interference_difference':float(np.sum(probs[mask])-np.sum(pinched[mask]))}
                 for name,mask in masks.items() if name!='bareP'}
    # Finite-window time average from the exact spectral sinc formula.
    de=prop.e[:,None]-prop.e[None,:]
    kernel=np.exp(1j*de*T/2)*np.sinc(de*T/(2*np.pi))
    cprod=prop.coeff.conj()[:,None]*prop.coeff[None,:]
    averages={}
    for name in ['q0','q2','R']:
        oe=prop.v.T@(masks[name][:,None]*prop.v)
        averages[name]=float(np.real(np.sum(cprod*oe*kernel)))
    low=prop.v[:,~high]@(phase[~high]*prop.coeff[~high]);low/=math.sqrt(1-w)
    normalized_low={name:float(np.sum(abs(low[mask])**2)) for name,mask in masks.items() if name!='bareP'}
    # Finite source/output supports of each low-energy eigenmode, not
    # labels of permanent trapping in an infinite-volume limit.
    low_modes=[]
    for i in np.flatnonzero(~high):
        low_modes.append({'energy':float(prop.e[i]),'prepared_weight':float(weights[i]),
            'both_fuel_mode_support':float(prop.v[prop.source,i]**2),
            'output_R_mode_support':float(np.sum(prop.v[masks['R'],i]**2)),
            'bareP_mode_support':float(np.sum(prop.v[masks['bareP'],i]**2)),
            'return_contribution_rel':pair(weights[i]*np.exp(-1j*(prop.e[i]-2*b1.OMEGA)*T))})
    return {'low_dimension':38,'high_energy_lower_bound':str(high_edge),
            'initial_mean_energy':str(2*Om),'initial_energy_variance':str(2*J*J),
            'high_band_weight_exact_upper':str(wbound),'high_band_weight_upper':float(wbound),
            'high_band_weight_computed':w,'probability_difference_from_normalized_low_upper':2*math.sqrt(2*float(wbound)),
            'last_low_eigenvalue':float(prop.e[37]),'first_high_eigenvalue':float(prop.e[38]),
            'eigenvalue_interval_radius_account':eta,'spectral_groups':len(groups),
            'low_spectral_groups':low_groups,'high_group_return_amplitude_rel':pair(high_return),
            'effective_group_participation':1/sum(x*x for x in group_weights),
            'endpoint_comparisons':comparisons,'finite_time_average_0_to_6':averages,
            'normalized_low_endpoint':normalized_low,'all_low_modes':low_modes,
            'scope':'All modes retained in evolution. Pinching is explanatory, not an implemented protocol or infinite-time replacement.'}


def interior_channel():
    model=b1.PathModel(3,sites=(1,),resource=2);colored=b2.Colored(model);n=model.dim
    hc=colored.h.toarray().reshape(n,4,n,4)
    aa=hc[:,1,:,1];bb=hc[:,2,:,1]
    S=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
    eye=np.eye(4);PT=(eye+S)/2;PS=(eye-S)/2
    defect=float(norm(colored.h.toarray()-np.kron(aa,eye)-np.kron(bb,S)))
    assert defect<1e-14 and norm(aa+bb-model.H.toarray())<1e-14
    src=model.index[((0,2),(0,2),0)]
    t=6/b1.JPAIR
    propT=b2.Propagator(csr_matrix(aa+bb),src);propS=b2.Propagator(csr_matrix(aa-bb),src)
    phiT=propT.state(t);phiS=propS.state(t);g=np.vdot(phiS,phiT)
    kraus=[a*PT+b*PS for a,b in zip(phiT,phiS)]
    tp=sum(k.conj().T@k for k in kraus)
    assert norm(tp-eye)<1e-11
    phase=np.exp(-1j*np.angle(g))*PT+PS
    fidelity=float(sum(abs(np.trace(phase@k))**2 for k in kraus)/16)
    assert abs(fidelity-(5+3*abs(g))/8)<1e-11
    def chan(rho):return sum(k@rho@k.conj().T for k in kraus)
    rho01=np.outer(eye[:,1],eye[:,1]);rho10=np.outer(eye[:,2],eye[:,2])
    td=float(np.sum(abs(np.linalg.eigvalsh(chan(rho01)-chan(rho10))))/2)
    assert abs(td-abs(g))<1e-11
    triplet=(eye[:,0]+1j*eye[:,3])/math.sqrt(2)
    assert norm(chan(np.outer(triplet,triplet.conj()))-np.outer(triplet,triplet.conj()))<1e-11
    pT=abs(phiT[src])**2;pS=abs(phiS[src])**2
    returned_average=(3*pT+pS)/4
    returned_fidelity=(3*math.sqrt(pT)+math.sqrt(pS))**2/(4*(3*pT+pS))
    herald=float(np.sum(np.minimum(abs(phiT)**2,abs(phiS)**2)))
    fuel=np.array([ff for p,m,ff in model.states])
    fuelblocks=[]
    for f in [0,1]:
        mask=fuel==f;pt=float(np.sum(abs(phiT[mask])**2));ps=float(np.sum(abs(phiS[mask])**2))
        gg=np.vdot(phiS[mask],phiT[mask])
        fuelblocks.append({'fuel':f,'probability_triplet':pt,'probability_singlet':ps,
                          'cross_gram':pair(gg),'unobserved_environment_cauchy_deficit':pt*ps-abs(gg)**2})
    rows=[]
    for i,(p,m,f) in enumerate(model.states):
        rows.append({'plus':p,'minus':m,'fuel':f,'amplitude_triplet':pair(phiT[i]),
                     'amplitude_singlet':pair(phiS[i]),'fields':model.fields[i]})
    return {'fixture':'conditional co-pairs at outer events, empty middle fuel; not original B2 fuel preparation',
            'theta':6,'time':t,'spatial_dimension':n,'colored_dimension':len(colored.states),
            'I_SWAP_representation_defect':defect,'gamma':pair(g),'gamma_modulus':float(abs(g)),
            'color_only_trace_distance_for_orthogonal_01_10':td,
            'entangled_input_fidelity_after_known_collective_phase':fidelity,
            'triplet_encoded_algebra':'exactly preserved by the analytic channel',
            'return_event_probability_triplet':float(pT),'return_event_probability_singlet':float(pS),
            'return_event_probability_maximally_mixed_input':float(returned_average),
            'return_event_conditional_entanglement_fidelity_after_phase':float(returned_fidelity),
            'full_configuration_measurement_and_filter_identity_success':herald,
            'fuel_record_gram_blocks':fuelblocks,
            'fuel_occupation_sector_distinguishability':abs(fuelblocks[1]['probability_triplet']-fuelblocks[1]['probability_singlet']),
            'TP_defect':float(norm(tp-eye)),'triplet_propagation_account':propT.account(),
            'singlet_propagation_account':propS.account(),'all_environment_amplitudes':rows,
            'scope':'Explicit diagnostic/recovery formulas, not an autonomous recovery or a replacement B2 success.'}


def main(out):
    begin=time.perf_counter();out.mkdir(parents=True,exist_ok=True)
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())['sources']
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in pins.items())
    model,prop,states,table,columns,sectors,masks=physical_account()
    spectrum=spectral_account(model,prop,states,masks);information=interior_channel()
    files=[]
    for name,data in [('SECTORS.json',sectors),('SPECTRUM.json',spectrum),('INFORMATION.json',information)]:
        (out/name).write_text(json.dumps(data,indent=2)+'\n');files.append(name)
    for a,start in enumerate(range(0,193,32)):
        name=f'RAW_{a}.json'
        (out/name).write_text(json.dumps({'columns':columns,'rows':table[start:start+32].tolist()},separators=(',',':'))+'\n')
        files.append(name)
    elapsed=time.perf_counter()-begin;rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    assert elapsed<60 and rss<1024
    print(json.dumps({'elapsed_seconds':elapsed,'peak_rss_MiB':rss,'data_files':files,
                      'data_sha256':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in files},
                      'historical_source_pins_unchanged':len(pins),'law':'unchanged B2/B1; no replacement test'},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out-dir',type=pathlib.Path,default=HERE)
    main(parser.parse_args().out_dir)

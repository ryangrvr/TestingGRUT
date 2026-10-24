#!/usr/bin/env python3
"""Conversion/redshift with an absorbing kinematic boundary and free receiver.

The momentum source is conditional on the parent's imposed scalar-at-rest
source and fixed ground objects. Path lengths are not a matter transfer fit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_core_conversion_20261023'
POINT=ROOT.parent/'grut_daughter_core_20261023/PARAMETER_POINT.json'
FIXED=ROOT.parent/'grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json'
H0=67.4*1000/3.085677581491367e22
OR=9.20917977964354e-5;OM=.315;OL=1-OM-OR
C=299792458./3.085677581491367e22  # Mpc/s
YQ=4.8294009076988175e-9;S0=2891.2;QOBJ=200000.
T0=2.7255*8.617333262e-14
GS0=2+7/8*6*(4/11)


def hubble(a):return H0*np.sqrt(OR/a**4+OM/a**3+OL)


def travel(a0,u0,af=1.,radiation=False):
    if a0>=af:return 0.
    q=a0*u0
    def integrand(x):
        a=np.exp(x);u=q/a
        h=H0*np.sqrt(OR)/a**2 if radiation else hubble(a)
        return C*u/np.sqrt(1+u*u)/(a*h)
    return quad(integrand,np.log(a0),np.log(af),epsabs=1e-13,epsrel=2e-9,limit=150)[0]


def map_weights(grid,locations):
    i=np.clip(np.searchsorted(grid,locations)-1,0,len(grid)-2)
    w=np.clip((locations-grid[i])/(grid[i+1]-grid[i]),0,1)
    return i,w


def mapped(values,i,w,n):
    return np.bincount(i,weights=values*(1-w),minlength=n)+np.bincount(i+1,weights=values*w,minlength=n)


def run(label,aform,high=400,low=150,ds_max=.001,source_width=.005,multiplier=1.,reactions=True):
    start=time.monotonic();packet=json.loads((PARENT/'RATE_RESULTS.json').read_text())
    omega=packet['Omega'];B=1+2*omega;uth=np.sqrt((B-1)*(B+1));e0=np.sqrt(1.5);u0=np.sqrt(.5)
    k0=e0-1
    grid=np.unique(np.r_[0.,np.geomspace(1e-18,1e-4,low),np.linspace(1e-4,k0+.002,high)])
    n=len(grid);kp=2*omega+grid;up=np.sqrt(kp*(kp+2));ua=np.sqrt(grid*(grid+2))
    ep=1+kp;ea=1+grid
    tab_u=np.array([r['p_over_m'] for r in packet['rows']]);tab_sv=np.array([r['sigma_v_conversion_cm3_s'] for r in packet['rows']])
    fit=PchipInterpolator(np.log(tab_u),np.log(tab_sv),extrapolate=False)
    sv=np.exp(fit(np.log(np.clip(ua,tab_u[0],tab_u[-1]))))*multiplier
    rev=sv*np.divide(ua*ea,up*ep,out=np.zeros_like(ea),where=up>0)
    if not reactions:sv[:]=0;rev[:]=0
    # Rows: positive/negative number, their first and second traveled-path
    # moments, and the population of initial anti-X which never converted.
    state=np.zeros((7,n));i,w=map_weights(grid,np.array([k0]));state[1,i[0]]=1-w[0];state[1,i[0]+1]=w[0];state[6]=state[1]
    nsrc=int(np.ceil(np.log(1/aform)/source_width))+2
    # Per source bin: N, N log(a_escape), prefix path first/second moments.
    source=np.zeros((4,nsrc));active_work=0.;net=0.;forward=0.;inverse=0.;history=[]
    maxN=0.;maxQ=0.;min_state=0.;snapshots=[]
    endpoints=[('equality',OR/OM),('z5',1/6),('today',1.)]
    nobj0=S0*YQ/QOBJ

    def half_drift(a,dh):
        nonlocal state,active_work,source
        # Stable kinetic energy maps retain tiny anti-X momenta explicitly.
        ua2=grid*(grid+2)*np.exp(-2*dh)
        newka=ua2/(np.sqrt(1+ua2)+1)
        newkp=(up*up*np.exp(-2*dh)-uth*uth)/(np.sqrt(1+up*up*np.exp(-2*dh))+B)
        escape=newkp<=0;stay=~escape
        # Integrate traveled path until a crossing, or the full half step.
        duration=np.minimum(dh,np.maximum(np.log(up/uth),0.))
        amidp=a*np.exp(duration/2);umidp=up*np.exp(-duration/2)
        dp=duration*C*umidp/np.sqrt(1+umidp*umidp)/(amidp*hubble(amidp))
        amida=a*np.exp(dh/2);umida=ua*np.exp(-dh/2)
        da=dh*C*umida/np.sqrt(1+umida*umida)/(amida*hubble(amida))
        state[4]+=2*dp*state[2]+dp*dp*state[0];state[2]+=dp*state[0]
        state[5]+=2*da*state[3]+da*da*state[1];state[3]+=da*state[1]
        active_work+=float(np.dot(state[0],kp-np.where(escape,2*omega,2*omega+newkp))
                           +np.dot(state[1],grid-newka))
        ae=a*np.exp(duration[escape]);weights=state[0,escape]
        index=np.clip(((np.log(ae/aform))/source_width).astype(int),0,nsrc-1)
        for row,val in enumerate([weights,weights*np.log(ae),state[2,escape],state[4,escape]]):
            source[row]+=np.bincount(index,weights=val,minlength=nsrc)
        ia,wa=map_weights(grid,newka);ip,wp=map_weights(grid,newkp[stay])
        for row in [0,2,4]:state[row]=mapped(state[row,stay],ip,wp,n)
        for row in [1,3,5,6]:state[row]=mapped(state[row],ia,wa,n)

    def snapshot(a,name):
        valid=(source[0]>1e-17)
        count=source[0,valid];ae=np.exp(source[1,valid]/count)
        pre1=source[2,valid];pre2=source[3,valid]
        lengths=np.array([travel(x,uth,a) for x in ae])
        q=uth*ae;u=q/a;energy=np.sqrt(1+u*u)
        mean1=float(np.dot(count,lengths));mean2=float(np.dot(count,lengths*lengths))
        total1=float(sum(pre1+count*lengths)+sum(state[2:4].ravel()))
        total2=float(sum(pre2+2*pre1*lengths+count*lengths*lengths)+sum(state[4:6].ravel()))
        original_length=travel(aform,u0,a)
        total1+=original_length;total2+=original_length**2
        original_energy=np.sqrt(1+(aform*u0/a)**2)
        source_energy=float(np.dot(count,energy))
        omitted=float(sum(source[0])-sum(count))
        # Omitted receiver weight is negligible and remains disclosed.
        receiver_work=float(np.dot(count,B-energy))
        Efree=float(original_energy+source_energy+np.dot(state[0],ep)+np.dot(state[1],ea))
        work=active_work+receiver_work+e0-original_energy
        residual=(Efree-2*omega*net+work-2*e0)/(2*e0)
        return dict(name=name,a=float(a),z=float(1/a-1),N_converted_receiver=float(sum(source[0])),
            omitted_source_weight=omitted,active_positive=float(sum(state[0])),active_anti=float(sum(state[1])),
            first_conversion_fraction=float(1-sum(state[6])),net_conversion=net,
            forward_events=forward,inverse_events=inverse,
            converted_post_escape_mean_Mpc=mean1/max(sum(count),1e-100),
            converted_post_escape_rms_Mpc=np.sqrt(mean2/max(sum(count),1e-100)),
            all_X_total_path_mean_Mpc=total1/2,all_X_total_path_rms_Mpc=np.sqrt(total2/2),
            original_X_collisionless_Mpc=original_length,
            receiver_mean_p_over_m=float(np.dot(count,u)/max(sum(count),1e-100)),
            receiver_number_fraction_p_over_m_above_01=float(sum(count[u>.01])),
            receiver_number_fraction_p_over_m_above_1=float(sum(count[u>1.])),
            converted_number_with_travel_above_Mpc={str(cut):float(sum(count[lengths>cut])) for cut in [.1,1.,5.,10.]},
            Efree_over_m_per_initial_pair=Efree,expansion_work_over_m=work,
            core_chemical_energy_change_over_m=-2*omega*net,
            energy_ledger_relative_error=float(residual),
            source_bins=[dict(a=float(x),N=float(w),q_over_m=float(uth*x),
                 post_escape_to_observation_Mpc=float(l),prefix_path_mean_Mpc=float(m/w))
                 for x,w,l,m in zip(ae,count,lengths,pre1)])

    current=aform
    for name,end in endpoints:
        if end<=current:continue
        segment=np.log(end/current);steps=max(1,int(np.ceil(segment/ds_max)));ds=segment/steps
        beginning=current
        for step in range(steps):
            a=beginning*np.exp(step*ds)
            half_drift(a,ds/2)
            amid=a*np.exp(ds/2);A=nobj0/amid**3*sv/hubble(amid);D=nobj0/amid**3*rev/hubble(amid)
            decay=A+D;f=np.divide(A,decay,out=np.zeros_like(A),where=decay>0);exp=np.exp(-decay*ds)
            oldplus=state[0].copy();total=state[0]+state[1];eq=f*total
            transient=np.divide(-np.expm1(-decay*ds),decay,out=np.full_like(decay,ds),where=decay>0)
            intplus=eq*ds+(oldplus-eq)*transient;intanti=total*ds-intplus
            forward+=float(np.dot(A,intanti));inverse+=float(np.dot(D,intplus))
            for p,m in [(0,1),(2,3),(4,5)]:
                total=state[p]+state[m];eq=f*total
                state[p]=eq+(state[p]-eq)*exp;state[m]=total-state[p]
            state[6]*=np.exp(-A*ds)
            net+=float(sum(state[0]-oldplus))
            half_drift(amid,ds/2)
            last=a*np.exp(ds)
            number=1+sum(source[0])+sum(state[:2].ravel())
            charge=1+sum(source[0])+sum(state[0])-sum(state[1])-2*net
            maxN=max(maxN,abs(number-2));maxQ=max(maxQ,abs(charge))
            min_state=min(min_state,float(np.min(state[[0,1,6]])))
            if step%max(1,int(.05/ds))==0 or step+1==steps:
                history.append([float(last),float(1-sum(state[6])),float(sum(state[1])),float(sum(source[0])),net,forward,inverse])
            if sum(state[:2].ravel())<1e-18:
                # No appreciable reactive population remains; the receiver
                # supplies exact later collisionless momenta and horizons.
                break
        current=end;snapshots.append(snapshot(end,name))
    # A final source CDF provides the last inelastic escape history.
    valid=source[0]>1e-17;count=source[0,valid];ae=np.exp(source[1,valid]/count)
    cdf=np.cumsum(count)/max(sum(count),1e-100)
    quantiles={str(q):float(1/np.interp(q,cdf,ae)-1) for q in [.1,.5,.9,.99]} if len(ae) else {}
    escape_source=snapshots[-1].pop('source_bins')
    for snap in snapshots:snap.pop('source_bins',None)
    return dict(scope='Conditional imposed source and fixed ground objects. Collisionless receiver horizons and total path lengths are not a transfer-function or Lyman-alpha likelihood.',
        label=label,a_form=aform,z_form=1/aform-1,omega_over_m=omega,
        escape_energy_over_m=B,escape_p_over_m=uth,nodes=n,ds_max=ds_max,source_bin_width=source_width,
        rate_multiplier=multiplier,reactions=reactions,first_conversion_history=history,
        last_inelastic_escape_redshift_quantiles=quantiles if sum(count)>0 else {},
        maximum_number_absolute_error=maxN,maximum_charge_absolute_error=maxQ,minimum_population=min_state,
        rate_table_sha256=hashlib.sha256((PARENT/'RATE_RESULTS.json').read_bytes()).hexdigest(),
        snapshots=snapshots,escape_source=escape_source,elapsed_seconds=time.monotonic()-start)


def mass_accounting():
    d=json.loads(FIXED.read_text());pop=d['rows'][1]['joint_screen_population'];point=json.loads(POINT.read_text())
    omega=json.loads((PARENT/'RATE_RESULTS.json').read_text())['Omega'];mx=point['m_X_GeV']
    pairs=pop['late_X_particles_per_object']/2;mass=d['rows'][1]['mass_GeV'];total=pop['late_mass_per_object_GeV']
    release=2*omega*mx*pairs;remaining=total-release
    return dict(m_X_GeV=mx,omega_over_m=omega,E_over_Qm=d['rows'][1]['energy_per_free_mass'],
        initial_anti_fraction_of_calibrated_matter=mx*pairs/total,
        complete_net_conversion_core_charge_fraction=2*pairs/QOBJ,
        complete_net_conversion_core_energy_GeV=release,
        released_energy_fraction_of_original_calibrated_rest_matter=release/total,
        remaining_rest_matter_ratio_without_retuning=remaining/total,
        final_converted_X_fraction_of_remaining_rest_matter=mx*pairs/remaining,
        final_all_free_X_fraction_of_remaining_rest_matter=2*mx*pairs/remaining,
        cold_ejected_K_over_m=2*omega,cold_ejected_p_over_m=np.sqrt((1+2*omega)**2-1),
        initial_daughter_ejected_K_over_m=np.sqrt(1.5)-1+2*omega,
        initial_daughter_ejected_p_over_m=np.sqrt((np.sqrt(1.5)+2*omega)**2-1),
        pairs_per_reference_object=pairs,reference_object_mass_GeV=mass,reference_late_total_mass_GeV=total)


def main():
    p=argparse.ArgumentParser();p.add_argument('--quick',action='store_true');args=p.parse_args()
    a20=T0/.02*(GS0/10.75)**(1/3);akev=T0/1e-6
    cases=[('T20MeV',a20),('T1keV',akev),('z1e6',1/1000001)]
    rows=[]
    settings=[(200,90,.002,.01)] if args.quick else [(200,90,.002,.01),(400,150,.001,.005),(800,250,.0005,.0025)]
    for label,af in cases:
        for high,low,ds,width in settings:
            r=run(label,af,high,low,ds,width);rows.append(r)
            last=r['snapshots'][-1];print(label,r['nodes'],'seconds',r['elapsed_seconds'],
              'Nrecv',last['N_converted_receiver'],'path Mpc',last['converted_post_escape_mean_Mpc'],flush=True)
            (ROOT/f"RUN_{label}_{r['nodes']}.json").write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
            compact=[{k:v for k,v in row.items() if k not in ['first_conversion_history','escape_source']} for row in rows]
            (ROOT/'STREAMING_RESULTS.json').write_text(json.dumps(dict(mass=mass_accounting(),runs=compact),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()

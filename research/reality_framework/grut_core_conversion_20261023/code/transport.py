#!/usr/bin/env python3
"""Conditional two-species conversion/redshift network, fixed cold objects.

This is a transport method test using the parent's imposed slow-s source.
It does not consume a finished formation source. The profile is frozen;
the core charge and its linear chemical energy are retained as ledgers.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scattering import ROOT,MX

H0=67.4*1000/3.085677581491367e22
OR=9.20917977964354e-5
OM=.315
YQ=4.8294009076988175e-9
S0=2891.2  # cm^-3; supplied present entropy normalization, held fixed
QOBJ=200000.


def rates(nu,omega,table,multiplier=1.):
    E=nu-omega;u=np.sqrt(np.maximum(E*E-1,0.))
    tab_u=np.array([r['p_over_m'] for r in table])
    tab_sv=np.array([r['sigma_v_conversion_cm3_s'] for r in table])
    fit=PchipInterpolator(np.log(tab_u),np.log(tab_sv),extrapolate=False)
    value=np.exp(fit(np.log(np.clip(u,tab_u[0],tab_u[-1]))))*multiplier
    value[E<1]=0.
    # Detailed balance: k_-^2 sigma_-+ = k_+^2 sigma_+-.
    Ep=nu+omega;up=np.sqrt(np.maximum(Ep*Ep-1,0.))
    reverse=value*np.divide(u*E,up*Ep,out=np.zeros_like(E),where=up>0)
    return value,reverse


def deposit(grid,values,location):
    i=np.clip(np.searchsorted(grid,location)-1,0,len(grid)-2)
    w=(location-grid[i])/(grid[i+1]-grid[i]);w=np.clip(w,0,1)
    return np.bincount(i,weights=values*(1-w),minlength=len(grid))+np.bincount(i+1,weights=values*w,minlength=len(grid))


def simulate(nodes=801,ds=.001,formation_z=1e6,multiplier=1.,packet=None):
    if packet is None:packet=json.loads((ROOT/'RATE_RESULTS.json').read_text())
    omega=packet['Omega']
    e0=np.sqrt(1.5);threshold=1+omega
    low=np.linspace(1-omega,threshold,nodes//2+1)
    high=np.linspace(threshold,e0+omega+.002,nodes//2+1)
    grid=np.r_[low,high[1:]]
    # One initial X and one initial anti-X; report per initial pair.
    def mono(location):
        out=np.zeros_like(grid);i=np.clip(np.searchsorted(grid,location)-1,0,len(grid)-2)
        w=(location-grid[i])/(grid[i+1]-grid[i]);out[i]=1-w;out[i+1]=w;return out
    plus=mono(e0-omega);anti=mono(e0+omega)
    sv,srev=rates(grid,omega,packet['rows'],multiplier)
    ap=1/(1+formation_z);span=np.log(1/ap);steps=int(np.ceil(span/ds));ds=span/steps
    # Exact free redshift maps over half a split step; interpolation preserves
    # particle number and mean energy of the represented points.
    Ep=grid+omega;Ea=np.maximum(grid-omega,1.)
    map_p=np.sqrt(1+(Ep*Ep-1)*np.exp(-ds))-omega
    map_a=np.sqrt(1+(Ea*Ea-1)*np.exp(-ds))+omega
    nobj0=S0*YQ/QOBJ
    net_conversion=0.;forward_events=0.;inverse_events=0.;redshift_loss=0.
    Ninitial=float(sum(plus+anti));Einitial=float(np.dot(plus,Ep)+np.dot(anti,Ea))
    maximum_N_error=0.;maximum_energy_ledger_error=0.;hist=[]
    for step in range(steps):
        a=ap*np.exp((step+.5)*ds)
        hubble=H0*np.sqrt(OR/a**4+OM/a**3+1-OM-OR)
        before=np.dot(plus,Ep)+np.dot(anti,Ea)
        plus=deposit(grid,plus,map_p);anti=deposit(grid,anti,map_a)
        after=np.dot(plus,Ep)+np.dot(anti,Ea);redshift_loss+=before-after
        A=nobj0/a**3*sv/hubble;B=nobj0/a**3*srev/hubble
        decay=A+B;total=plus+anti
        eq=np.divide(A*total,decay,out=plus.copy(),where=decay>0)
        change=(eq-plus)*(-np.expm1(-decay*ds))
        # Integrals of each state through the exact 2x2 reaction update.
        transient=np.divide(-np.expm1(-decay*ds),decay,out=np.full_like(decay,ds),where=decay>0)
        intplus=eq*ds+(plus-eq)*transient;intanti=total*ds-intplus
        forward_events+=float(np.dot(A,intanti));inverse_events+=float(np.dot(B,intplus))
        plus+=change;anti-=change;net_conversion+=float(sum(change))
        before=np.dot(plus,Ep)+np.dot(anti,Ea)
        plus=deposit(grid,plus,map_p);anti=deposit(grid,anti,map_a)
        after=np.dot(plus,Ep)+np.dot(anti,Ea);redshift_loss+=before-after
        maximum_N_error=max(maximum_N_error,abs(sum(plus+anti)/Ninitial-1))
        residual=after-2*omega*net_conversion+redshift_loss-Einitial
        maximum_energy_ledger_error=max(maximum_energy_ledger_error,abs(residual)/Einitial)
        if step%max(1,steps//200)==0 or step+1==steps:
            hist.append([float(ap*np.exp((step+1)*ds)),float(sum(plus)),float(sum(anti)),float(net_conversion),
                         float(redshift_loss),float(after)])
    return dict(scope='Imposed monoenergetic slow-s daughters, identical fixed ground objects, linear conversion and redshift only; no actual formation, mergers, nonlinear capture or completed cosmology.',
        rate_table_numeric_sha256=hashlib.sha256(json.dumps(packet['rows'],sort_keys=True).encode()).hexdigest(),
        nodes=len(grid),ds=ds,steps=steps,formation_z=formation_z,rate_multiplier=multiplier,
        remaining_anti_fraction_of_initial_anti=float(sum(anti)),net_anti_conversion=net_conversion,
        forward_events_per_initial_pair=forward_events,inverse_events_per_initial_pair=inverse_events,
        number_relative_error=maximum_N_error,energy_ledger_relative_error=maximum_energy_ledger_error,
        free_charge_increase_per_initial_pair=float(sum(plus-anti)),
        core_charge_change_per_initial_pair=-2*net_conversion,
        core_chemical_energy_change_over_mX_per_initial_pair=-2*omega*net_conversion,
        source_energy_over_mX_per_initial_pair=Einitial,
        final_free_energy_over_mX_per_initial_pair=float(np.dot(plus,Ep)+np.dot(anti,Ea)),
        redshift_work_over_mX_per_initial_pair=redshift_loss,
        history=hist)


def main():
    rows=[];packet=json.loads((ROOT/'RATE_RESULTS.json').read_text())
    for z in [1e6,1e8,1.2e11]:
        for n,ds in [(401,.002),(801,.001),(1601,.0005)]:
            row=simulate(n,ds,z,packet=packet);rows.append(row)
            print(z,n,'anti survival',row['remaining_anti_fraction_of_initial_anti'],
                  'forward/inverse',row['forward_events_per_initial_pair'],row['inverse_events_per_initial_pair'],flush=True)
    for mult in [.1,10.]:
        row=simulate(801,.001,1.2e11,mult,packet=packet);rows.append(row)
    (ROOT/'TRANSPORT_RESULTS.json').write_text(json.dumps(dict(rows=rows),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()

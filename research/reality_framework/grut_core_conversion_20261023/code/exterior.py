#!/usr/bin/env python3
"""Measured flux ratios and conditional equal-time exterior mode projections.

The masked snapshot is not a complete asymptotic escaping source. The
absorbed field and unresolved low-k band are never filled in by normalization.
"""
import json
from pathlib import Path
import numpy as np
from scipy.fft import dst
from scattering import ROOT

PARENT=ROOT.parent/'grut_daughter_core_20261023'


def taper(r,start,width):
    t=np.clip((r-start)/width,0,1)
    return np.sin(np.pi*t/2)**2


def project(table,parameters,start=80.,width=20.):
    r=table[:,0];h=parameters['h'];L=parameters['length'];nodal=parameters.get('nodal',False)
    mask=taper(r,start,width)*(1-taper(r,parameters['absorber']-40,40))
    psi=table[:,1]+1j*table[:,2];pi=table[:,4]+1j*table[:,5]
    u=mask*r*psi;p=mask*r*pi
    kind=1 if nodal else 2
    U=dst(u.real,type=kind,norm='ortho')+1j*dst(u.imag,type=kind,norm='ortho')
    P=dst(p.real,type=kind,norm='ortho')+1j*dst(p.imag,type=kind,norm='ortho')
    k=np.pi*np.arange(1,len(r)+1)/L;omega=np.sqrt(1+k*k)
    plus=4*np.pi*h*abs((np.sqrt(omega)*U+1j*P/np.sqrt(omega))/2)**2
    minus=4*np.pi*h*abs((np.sqrt(omega)*U-1j*P/np.sqrt(omega))/2)**2
    energy=omega*(plus+minus)
    active=mask>.01
    correction=abs(table[:,3]**2-1+.1*abs(psi)**2)
    change=abs(2*table[:,3]*table[:,6]+.2*np.real(psi.conj()*pi))
    delta=float(np.max(correction[active]));deltadot=float(np.max(change[active]))
    controlled=k*k>max(10*delta,1e-4)
    qspace=-4*np.pi*h*np.sum(np.imag(u.conj()*p))
    return dict(inner_taper_start=start,inner_taper_width=width,outer_taper_end=parameters['absorber'],
        max_local_mass2_correction=delta,max_time_derivative_mass2=deltadot,
        Nplus_bar=float(sum(plus)),Nminus_bar=float(sum(minus)),Qbar=float(sum(plus-minus)),
        charge_parseval_residual=float((sum(plus-minus)-qspace)/max(abs(qspace),1e-100)),
        Ebar_free=float(sum(energy)),mQ_over_E=float(qspace/sum(energy)),
        conditional_antiparticle_number_fraction=float(sum(minus)/sum(plus+minus)),
        controlled_free_energy_fraction=float(sum(energy[controlled])/sum(energy)),
        controlled_Nplus_bar=float(sum(plus[controlled])),controlled_Nminus_bar=float(sum(minus[controlled])),
        unresolved_free_energy_fraction=float(sum(energy[~controlled])/sum(energy)))


def main():
    spectra=json.loads((PARENT/'SPECTRAL_RESULTS.json').read_text())['rows']
    flux=[]
    for row in spectra:
        x=row['X'];E=x['total_windowed_net_energy'];Q=x['total_windowed_net_charge'];ratio=Q/E
        flux.append(dict(run=row['run'],Ebar_net_windowed=E,Qbar_net_windowed=Q,mQ_over_E=ratio,
            conditional_upper_antiparticle_number_fraction=(1-ratio)/2,
            controlled_Nplus_bar=x['controlled_positive_frequency_Nbar'],
            controlled_Nminus_bar=x['controlled_negative_frequency_Nbar'],
            declared_controlled_energy_fraction=x['controlled_fraction_of_outgoing_linear_energy'],
            incoming_over_outgoing_linear_energy=x['incoming_linear_energy']/x['outgoing_linear_energy']))
    projections=[]
    runs=json.loads((PARENT/'CORE_RESULTS.json').read_text())['runs']
    for run in runs:
        table=np.loadtxt(PARENT/run['final_file'],delimiter=',',skiprows=1)
        projections.append(dict(run='parent_'+run['name'],time=run['parameters']['tmax'],
            variants=[project(table,run['parameters'],s,w) for s,w in [(60,20),(80,20),(80,40),(100,20)]]))
    for path in sorted(ROOT.glob('*_RESULTS.json')):
        if not path.name.startswith(('H','N')):continue
        run=json.loads(path.read_text());table_path=ROOT/'data'/f"{run['name']}_final.csv"
        if not table_path.exists():continue
        table=np.loadtxt(table_path,delimiter=',',skiprows=1)
        projections.append(dict(run=run['name'],time=run['parameters']['tmax'],
            variants=[project(table,run['parameters'],s,w) for s,w in [(60,20),(80,20),(80,40),(100,20)]]))
    out=dict(scope='Flux ratios are measured. Snapshot particle numbers assume a free exterior and depend on masks; no complete asymptotic source or recovery of absorbed particles is claimed.',
        flux_ratio_caution='Kinetic energy alone lowers mQ/E. The antiparticle upper bound requires a nonnegative free outgoing ensemble, not arbitrary interacting net flux.',
        flux_rows=flux,exterior_snapshots=projections)
    (ROOT/'EXTERIOR_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for row in projections:
        v=row['variants'][1]
        print(row['run'],'snapshot E',v['Ebar_free'],'controlled fraction',v['controlled_free_energy_fraction'],
              'anti fraction',v['conditional_antiparticle_number_fraction'])


if __name__=='__main__':main()

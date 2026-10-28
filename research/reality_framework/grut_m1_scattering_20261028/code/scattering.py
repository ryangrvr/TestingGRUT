"""Flux-unitary three-channel radial M1 scattering, on a refined FD2 grid.

Variables are r times (u,v,sqrt(2)*w), where delta z = u exp(-i nu t)
+v* exp(+i nu t), delta S=w exp(-i nu t)+w* exp(+i nu t), and z rotates
as exp(-i Omega t). Thus the radial potential matrix is real symmetric.
"""
from common import *
from scipy.linalg import solve_banded

def exterior_ratio(k2,ell,h,R,Rfar=500.):
    """Outgoing/decaying discrete free radial solution by a Riccati recurrence.

    q_i=y_(i+1)/y_i. Backward propagation carries the centrifugal tail into
    the matching boundary. For open channels Im(q)/h is exactly the discrete
    radial flux per boundary amplitude squared; evanescent q is real.
    """
    N=round(R/h);far=max(N+1,round(Rfar/h))
    if k2>0:
        if h*np.sqrt(k2)>=2:raise ValueError('Propagating channel beyond lattice band')
        phase=2*np.arcsin(.5*h*np.sqrt(k2));q=complex(np.exp(1j*phase))
    else:q=complex(np.exp(-2*np.arcsinh(.5*h*np.sqrt(-k2))))
    for i in range(far,N,-1):q=1/(2+h*h*(ell*(ell+1)/(i*h)**2-k2)-q)
    return q

def scatter(k,ell=0,h=.05,R=60.,Rfar=500.,control=None,return_field=False,incoming_index=1,return_amplitudes=False):
    if k<=0:raise ValueError('Use a positive incident momentum')
    N=round(R/h);r=np.arange(1,N+1)*h;F,S=profile(r)
    if control=='vacuum':F=np.zeros_like(F);S=np.ones_like(S)
    E=np.sqrt(1+k*k);nu=E+OMEGA
    energies=np.array([E+2*OMEGA,E,nu]);masses2=np.array([1.,1.,6.]);k2=energies**2-masses2
    q=np.array([exterior_ratio(v,ell,h,R,Rfar) for v in k2])
    opened=k2>0;flux=np.where(opened,q.imag/h,0.)
    if np.any(flux[opened]<=0):raise RuntimeError('Outgoing flux sign')
    normal=S*S+2*A*F*F
    neutral=B*(3*S*S-1)-6*C*(S-1)+F*F
    centrifugal=ell*(ell+1)/(r*r)
    diag=np.c_[normal-energies[0]**2,normal-energies[1]**2,neutral-nu*nu]+centrifugal[:,None]+2/h**2
    ab=np.zeros((7,3*N),dtype=complex)
    idx=np.arange(3*N);ab[3]=diag.ravel()
    # -1/h² between adjacent radial blocks, for the same channel.
    ab[0,3:]=-1/h**2;ab[6,:-3]=-1/h**2
    mixing=A*F*F;coupling=np.sqrt(2)*F*S
    if control in ('vacuum','uncoupled_well'):mixing=np.zeros_like(F);coupling=np.zeros_like(F)
    # ab[upper+i-j,j]=M[i,j]. Fill each real-symmetric within-block pair.
    for left,right,value in [(0,1,mixing),(0,2,coupling),(1,2,coupling)]:
        ii=np.arange(N)*3+left;jj=np.arange(N)*3+right
        ab[3+left-right,jj]=value;ab[3+right-left,ii]=value
    ab[3,-3:]-=q/h**2
    if incoming_index not in (0,1,2) or not opened[incoming_index]:raise ValueError('Incident channel must propagate')
    incoming=1/np.sqrt(flux[incoming_index]);rhs=np.zeros(3*N,dtype=complex)
    rhs[-3+incoming_index]=(np.conj(q[incoming_index])-q[incoming_index])*incoming/h**2
    solution=solve_banded((3,3),ab,rhs,check_finite=False).reshape(N,3)
    outgoing=solution[-1].copy();outgoing[incoming_index]-=incoming
    amplitudes=np.zeros(3,dtype=complex);amplitudes[opened]=outgoing[opened]*np.sqrt(flux[opened])
    probabilities=abs(amplitudes)**2
    result={'k':float(k),'ell':ell,'h':h,'R':R,'Rfar':Rfar,'incoming_energy':float(E),
            'outgoing_X_energy':float(energies[0]),'outgoing_phi_energy':float(nu),
            'P_X':float(probabilities[0]),'P_antiparticle':float(probabilities[1]),'P_phi':float(probabilities[2]),
            'unitarity_error':float(abs(probabilities.sum()-1)),
            'core_charge_change_per_incident':float(-2*probabilities[0]-probabilities[2]),
            'anti_survival_after_one_vacuum_phi_decay':float(probabilities[1]+probabilities[2]),
            'outgoing_energy_gain':float(np.dot(energies,probabilities)-E),
            'rotating_energy_balance':float(abs(np.dot(energies,probabilities)-E-OMEGA*(2*probabilities[0]+probabilities[2]))),
            'control':control}
    if incoming_index!=1:
        qin=[1.,-1.,0.][incoming_index]
        shift=qin-(probabilities[0]-probabilities[1])
        gain=float(np.dot(energies,probabilities)-energies[incoming_index])
        result.update(incoming_energy=float(energies[incoming_index]),core_charge_change_per_incident=float(shift),
                      outgoing_energy_gain=gain,rotating_energy_balance=abs(gain+OMEGA*shift),
                      anti_survival_after_one_vacuum_phi_decay=None)
    if return_amplitudes:result['amplitudes_real_imag']=[[float(x.real),float(x.imag)] for x in amplitudes]
    result['incident_channel']=('X','antiparticle','phi')[incoming_index]
    if return_field:return result,r,solution
    return result

def main():
    ks=[.2,.4,.7,1.,1.15,1.25,1.5,2.,3.];rows=[]
    for k in ks:
        row=scatter(k);rows.append(row);print(json.dumps(row),flush=True)
    controls=[scatter(1.,control=x) for x in ('vacuum','uncoupled_well')]
    refinement=[]
    for k in (.4,1.,1.25,2.):
        refinement.append({'k':k,'coarse':scatter(k,h=.1),'fine':scatter(k,h=.025),
                           'far_boundary':scatter(k,R=80.),'far_matching':scatter(k,Rfar=1000.)})
    save('SCATTERING_RESULTS.json',{'scope':__doc__,'rows':rows,'controls':controls,'refinements':refinement,
                                   'full_cosmic_capture_rate':None})

if __name__=='__main__':main()

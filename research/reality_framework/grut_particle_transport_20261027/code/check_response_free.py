"""The matrix-response convention on an exact free FD4 lattice oscillator."""
from pt_common import *

def main():
    N,L=8,5.;set_box(L);lib=library();z=np.zeros((6,N,N,N));d=np.zeros((3,6,N,N,N))
    source=np.broadcast_to(np.cos(2*np.pi*np.arange(N)/N)[:,None,None],(N,N,N)).copy()
    for a in range(3):d[a,3+a]=source
    _,kl=geometry(N,L);omega=np.sqrt(kl[1,0,0]);times=np.arange(151)*.04;R=[];previous=0.
    for t in times:
        if t:step(lib,z,t-previous,.004,tangent=d,threads=2)
        previous=t;R.append(np.einsum('baix,ix->ba',d[:,:3].reshape(3,3,N,-1),source.reshape(N,-1))/np.sum(source*source))
    R=np.array(R);expected=np.sin(omega*times)/omega
    error=float(max(np.max(abs(R[:,a,a]-expected)) for a in range(3)))
    off=R.copy()
    for a in range(3):off[:,a,a]=0
    result=dict(passed=error<2e-8,exact_frequency=omega,maximum_response_absolute_error=error,
        maximum_cross_response=float(np.max(abs(off))),maximum_anomalous_response=float(np.max(abs((R[:,0,0]-R[:,1,1])/2))),
        scope='Exact free-background test of the discrete tangent and normal/conjugate response convention; no particle assignment in the nonlinear state.')
    save('RESPONSE_FREE_CHECK.json',result);print(json.dumps(result,indent=2));assert result['passed']

if __name__=='__main__':main()

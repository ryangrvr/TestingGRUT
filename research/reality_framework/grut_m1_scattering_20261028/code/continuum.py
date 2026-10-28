"""Separate continuum collocation check of the monopole scattering problem."""
from common import *
from scipy.integrate import solve_bvp
from scattering import scatter

def solve(k,tol=2e-8,R=60.):
    E=np.sqrt(1+k*k);nu=E+OMEGA;energies=np.array([E+2*OMEGA,E,nu]);k2=energies**2-[1,1,6]
    momenta=np.sqrt(k2.astype(complex));ain=1/np.sqrt(k)
    def crhs(r,y):
        f,s=profile(r);normal=s*s+2*A*f*f;ds=B*(3*s*s-1)-6*C*(s-1)+f*f
        d=np.array([normal-energies[0]**2,normal-energies[1]**2,ds-nu*nu])
        acc=d*y[:3];acc[0]+=A*f*f*y[1]+np.sqrt(2)*f*s*y[2]
        acc[1]+=A*f*f*y[0]+np.sqrt(2)*f*s*y[2]
        acc[2]+=np.sqrt(2)*f*s*(y[0]+y[1])
        return np.r_[y[3:],acc]
    src=np.array([0,-2j*k*ain,0])
    def cbc(l,r):return np.r_[l[:3],r[3:]-1j*momenta*r[:3]-src]
    def rhs(r,y):
        z=crhs(r,y[:6]+1j*y[6:]);return np.r_[z.real,z.imag]
    def bc(l,r):
        z=cbc(l[:6]+1j*l[6:],r[:6]+1j*r[6:]);return np.r_[z.real,z.imag]
    x=np.linspace(0,R,1201);y=np.zeros((12,len(x)))
    sol=solve_bvp(rhs,bc,x,y,tol=tol,bc_tol=tol,max_nodes=25000)
    if not sol.success:raise RuntimeError(sol.message)
    outgoing=(sol.y[:3,-1]+1j*sol.y[6:9,-1]);outgoing[1]-=ain
    prob=abs(outgoing)**2*np.where(k2>0,momenta.real,0)
    return {'k':k,'P_X':float(prob[0]),'P_antiparticle':float(prob[1]),'P_phi':float(prob[2]),
            'unitarity_error':float(abs(prob.sum()-1)),'max_residual':float(sol.rms_residuals.max()),'nodes':len(sol.x)}

def main():
    out=[]
    for k in (.4,1.,1.25,1.5,2.):
        c=solve(k);fd=scatter(k,h=.025)
        out.append({'continuum':c,'FD025':fd,'P_X_difference':fd['P_X']-c['P_X'],
                    'P_phi_difference':fd['P_phi']-c['P_phi']});print(json.dumps(out[-1]),flush=True)
    save('CONTINUUM_CHECKS.json',{'rows':out,'same_author_separate_method':True})

if __name__=='__main__':main()

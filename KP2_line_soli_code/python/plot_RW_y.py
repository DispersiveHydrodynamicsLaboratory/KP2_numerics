import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Code for plotting the Whitham Equation Solution to
# the KP2 equation with initial conditions of
# a line soliton with a step in its
# amplitude and angle of propagation in the y-direction only

# TODO: rescale solution back to Kodama-land

def main():
    # Define soliton parameters
    # Following Parameters match Kodama, case (h)
    ap = 2*6**(2/3)       # Upper half plane (UHP) soliton's amplitude
    qp  =  1*(3/np.sqrt(2))**(1/3) # UHP soliton's angle of propagation
    am = 2*6**(2/3)       # LHP soliton's amplitude
    qm  = -1*(3/np.sqrt(2))**(1/3) # LHP soliton's angle of propagation
    sap = np.sqrt(ap)
    sam = np.sqrt(am)
    
    
    # Define (x,y,t)
    Lx = 40
    Ly = 25
    Nx = 2**7
    Ny = 2**7
    x  = np.linspace(-Lx,Lx,Nx)
    y  = np.linspace(-Ly,Ly,Ny)
    [X,Y] = np.meshgrid(x,y)
    t  = 6*np.sqrt(3/8)
    
    # Run check for consistencies
    # Will return 0 and readout if inconsistencies
    # Otherwise will return 1(2) for Case 1(2)
    # cn = condition_check(sal,ql,sar,qr,m);
    # if cn==0:
    #     return False
    
    # Find solution in terms of R1 and R2
    ans = RI_soln(sam,qm,sap,qp,X,Y,t);
    if ans:
        (R1,R2) = ans
    else:
        return False
    
    sa = (R1+R2)/2
    q  = (R1-R2)/2
    U  = soli(X,Y,t,sa**2,q,0)
    # Plot solution (soliton)
    plt.figure(1)
    plt.contourf(X,Y,U)
    plt.colorbar()
    plt.show()
    plt.savefig('soli.png')
    # # Plot solution (amplitude)
    # plt.figure(1)
    # plt.contourf(X,Y,sa**2)
    # plt.colorbar()
    # plt.show()
    # plt.savefig('ampl.png')
    # # Plot solution (angle)
    # plt.figure(2)
    # plt.contourf(X,Y,q)
    # plt.colorbar()
    # plt.show()
    # plt.savefig('angle.png')
    
    
    return True
    
def condition_check(sal,ql,sar,qr,m):
    # Checks to ensure an all-RW solution will be given
    # Initialize casenumber, verbosity
    casenum = 0;
    verbose = 1;
    
    # Define Riemann Invariants
    (R1l,R2l) = RIs(sal,ql)
    (R1r,R2r) = RIs(sar,qr)
    
    # Check based on m<0
    if m>0:
        print('Current cases only allow for m<0. Please consider flipping your axes.')
        return False
    
    # Check based on middle speeds (i.e. no RW interaction)
    # Actual checks (MM: might be able to reduce these?)
    # ms0 = ( (R1r==2/m+R2l) & (R2l==-1/m) ) # should be in other cases
    ms1 = ( (R1r+R2l>=0) & (R2l>=-1/m) & (R1r<=2/m+R2l) )
    ms2 = ( (R1r+R2l<=0) & (R2l<=-1/m) & (R1r>=2/m+R2l)  )
    # if verbose:
    #     print(' ms1: ', ms1,' ms2: ',ms2)
    #     print('R1r: ',R1r,' R2l: ',R2l,' -1/m: ',-1/m)
    if not ( ms1 or ms2 ):
        print('Middle speed violation.')
        print(' ms1: ', ms1,' ms2: ',ms2)
        print('R1r: ',R1r,' R2l: ',R2l,' -1/m: ',-1/m)
        return False
    # See which case we're in; check speed conditions
    qc = (ql+qr) >= (-1/m);
    if qc:
        if verbose:
            print('Case 1 expected')
            print('R1l: ',R1l,' R1r: ',R1r)
            print('R2l: ',R2l,' R2r: ',R2r)
        R1c = R1l<=R1r;
        R2c = R2l>=R2r;
        if ~(R1c & R2c):
            print('R.I. violation for Case 1');
            print('R1c (R1l<=R1r): ', R1c, ' R2c (R2l>=R2r): ', R2c);
            return False
        if verbose:
            print('Case 1 confirmed')
        casenum = 1;
        return casenum
    else:
        if verbose:
            print('Case 2 expected')
            print('R1l: ',R1l,' R1r: ',R1r)
            print('R2l: ',R2l,' R2r: ',R2r)
        R1c = R1l>R1r;
        R2c = R2l<R2r;
        if ~(R1c & R2c):
            print('R.I. violation for Case 2');
            print('R1c (R1l>R1r): ', R1c, ' R2c (R2l<R2r): ', R2c);
            return False
        if verbose:
            print('Case 2 confirmed')
        casenum = 2;
        return casenum
    
def RI_soln(sam,qm,sap,qp,X,Y,t):
    # Define Riemann Invariants on LHS and RHS of jump
    (R1m,R2m) = RIs(sam,qm)
    (R1p,R2p) = RIs(sap,qp)
    
    # Define characteristic speeds
    s1  = V1(R1m,R2m)
    s2  = V1(R1p,R2m)
    s3  = V2(R1p,R2m)
    s4  = V2(R1p,R2p)
    
    # Check to make sure speeds are ordered
    sc = (s1<=s2)&(s2<=s3)&(s3<=s4)
    if not sc:
        print('Error! Speeds are not ordered correctly')
        print('sc: ',sc)
        print('1->2: ',(s1<=s2),' 2->3: ',(s2<=s3),' 3->4: ',(s3<=s4))
        print('s1: ','{:.3f}'.format(s1),' s2: ','{:.3f}'.format(s2),' s3: ','{:.3f}'.format(s3),' s4: ','{:.3f}'.format(s4))
        return False
    
    # # Determine solution, based on case
    # if cn==1:
    #     print(cn)
    #     # f = f_+ here
    #     # g = g_- here
    # elif cn==2:
    #     print(cn)
    #     # f = f_- here
    #     # g = g_+ here
    # else:
    #     print('Error! Not a possible case.')
    #     return(False)
    
    
    
    R1 = R1m * (Y <= s1*t) +  f(Y/t,R2m) * (s1*t < Y) * (Y < s2*t) + R1p * (Y >= s2*t)
    
    R2 = R2m * (Y <= s3*t) + g(Y/t,R1p) * (s3*t < Y) * (Y < s4*t) + R2p * (Y >= s4*t)
    # R1 = 0; R2 = 0;
    return (R1,R2)
    
def RIs(sa,q):
    # Shorthand for determining Riemann Invariants
    # Designed to intake parameters sa = \sqrt{a} and q
    R1 = sa + q
    R2 = sa - q
    return (R1,R2)
    
def V1(R1,R2):
    # Characteristic speed for first RI
    V1 = 2/3*R1 - 4/3*R2
    return V1
    
def V2(R1,R2):
    # Characteristic speed for second RI
    V2 = 4/3*R1 - 2/3*R2
    return V2
    
def f(xi,R2m):
    # 1-Wave solution
    RW1 = 2*R2m+3*xi/2
    return RW1
    
def g(xi,R1p):
    # 2-Wave solution
    RW2 = 2*R1p-3*xi/2
    return RW2
    
def soli(X,Y,t,A,Q,ubar):
    THETA = X + Q*Y - (ubar + A/3 + Q**2)*t
    U = ubar + A*(1/np.cosh(np.sqrt(A/12)*THETA)**2)
    return U
    
if __name__=='__main__':
    sys.exit(main())
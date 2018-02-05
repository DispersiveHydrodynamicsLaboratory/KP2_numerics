import sys
import numpy as np
import scipy as sp

# Code for plotting the Whitham Equation Solution to
# the KP2 equation with initial conditions of
# a line soliton with a step in its
# amplitude and angle of propagation

def main():
    # Define soliton parameters
    sal = 1.75          # Left hand side (LHS) soliton's amplitude
    ql  = 2+sal       # LHS soliton's angle of propagation
    sar = 0.25          # RHS soliton's amplitude
    qr  = -1-sar      # LHS soliton's angle of propagation
    m   = -1         # Slope of line over which jump in parameters occurs
    
    # Define (x,y,t)
    Lx = 30
    Ly = 30
    Nx = 2**6
    Ny = 2**6
    x  = np.linspace(-Lx,Lx,Nx)
    y  = np.linspace(-Ly,Ly,Ny)
    t  = 5
    
    # # Put soliton amplitudes in usable form
    # sal = np.sqrt(al)
    # sar = np.sqrt(ar)
    
    # Run check for consistencies
    # Will return 0 and readout if inconsistencies
    # Otherwise will return 1(2) for Case 1(2)
    cn = condition_check(sal,ql,sar,qr,m);
    if cn==0:
        return False
    
    # Find solution in terms of R1 and R2
    ans = RI_soln(sal,ql,sar,qr,m,cn);
    if ans:
        (R1,R2) = ans
    else:
        return False
    
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
    
def RI_soln(sal,ql,sar,qr,m,cn):
    # Define Riemann Invariants on LHS and RHS of jump
    (R1l,R2l) = RIs(sal,ql)
    (R1r,R2r) = RIs(sar,qr)
    
    # Define characteristic speeds
    mfac = -m/(3*np.sqrt(1+m**2)) # Factor in front of all speeds
    W1l  = mfac*(R1l+2/m)*(R1l-2*R2l)
    W1r  = mfac*(R1r+2/m)*(R1r-2*R2l)
    W2l  = mfac*(R2l-2/m)*(R2l-2*R1r)
    W2r  = mfac*(R2r-2/m)*(R2r-2*R1r)
    
    # Check to make sure speeds are ordered
    Wc = (W1l<=W1r)&(W1r<=W2l)&(W2l<=W2r)
    if ~Wc:
        print('Error! Speeds are not ordered correctly')
        print('1->2: ',(W1l<=W1r),' 2->3: ',(W1r<=W2l),' 3->4: ',(W2l<=W2r))
        print('W1l: ','{:.3f}'.format(W1l),' W1r: ','{:.3f}'.format(W1r),' W2l: ','{:.3f}'.format(W2l),' W2r: ','{:.3f}'.format(W2r))
        return False
    
    # Determine solution, based on case
    if cn==1:
        print(cn)
        # f = f_+ here
        # g = g_- here
    elif cn==2:
        print(cn)
        # f = f_- here
        # g = g_+ here
    else:
        print('Error! Not a possible case.')
        return(False)
    # R1 = R1l * (alpha <= W1l*t) + ...
    #      f(alpha/t) * (W1l*t < alpha) * (alpha < W1r*t) + ...
    #      R1r * (alpha >= W1r*r)
    
    # R2 = R2l * (alpha <= W2l*t) + ...
    #      g(alpha/t) * (W2l*t < alpha) * (alpha < W2r*t) + ...
    #      R2r * (alpha >= W2r*r)
    R1 = 0; R2 = 0;
    return (R1,R2)
    
def RIs(sa,q):
    # Shorthand for determining Riemann Invariants
    # Designed to intake parameters sa = \sqrt{a} and q
    R1 = sa + q
    R2 = sa - q
    return (R1,R2)
    
if __name__=='__main__':
    sys.exit(main())
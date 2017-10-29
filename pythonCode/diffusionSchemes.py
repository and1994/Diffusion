# Numerical schemes for simulating diffusion for outer code diffusion.py

from __future__ import absolute_import, division, print_function
import numpy as np

# The linear algebra package for BTCS (for solving the matrix equation)
import scipy.linalg as la

def FTCS(phiOld, d, nt):
    
    """
    Diffusion of profile in phiOld using FTCS using non-dimensional diffusion
    coefficient, d.
    Spatial endpoints obtained through FTFC and FTBS as FTCS could not be used
    there. Furthermore, boundary conditions (gradient of phi = 0 at the
    boundary points) are imposed.
    """
    
    # arguments test
    if nt<=0:
        raise ValueError('Error in FTCS: Argument nt to FTCS should be > 0')
    if not(int(nt) == nt):
        raise ValueError('Error in FTCS:\
                         Argument nt to FTCS should be an integer')
    if not(isinstance(float(d),float) and float(d) > 0):
        raise TypeError('Error in FTCS:\
                        Argument d to FTCS should be a positive float')
    if not(isinstance(phiOld,np.ndarray)):
        raise TypeError('Error in FTCS:\
                        Argument phiOld to FTCS should be an array')
    
    nx = len(phiOld)
    
    # new time-step array for phi
    phi = phiOld.copy()
    
    # FTCS for all time steps
    for it in range(int(nt)):
        
        #initial endpoint for spatial coordinates (imposing boundary conditions)
        phi[0] = phiOld[0] + d*(phiOld[0] - 2*phiOld[1] + phiOld[2])
        phi[1]=phi[0]
        
        # spatial points, endpoints excluded
        for j in range(2,nx-2):
            phi[j] = phiOld[j] + d*(phiOld[j+1] - 2*phiOld[j] + phiOld[j-1])
            
       
        # final endpoint for spatial coordinates (imposing boundary conditions)
        phi[nx-1] = phiOld[0] + d*(phiOld[nx-3] - 2*phiOld[nx-2] + phiOld[nx-1])
        phi[nx-2] = phi[nx-1]
        # output to phiOld for the next time-step
        phiOld = phi.copy()
        
    return phi
    
def BTCS(phi, d, nt):
    
    """
    Diffusion of profile in phi using BTCS using non-dimensional diffusion
    coefficient d, assuming fixed value boundary conditions.
    """
    # arguments test
    if nt<=0:
        raise ValueError('Error in BTCS: Argument nt to BTCS should be > 0')
    if not(int(nt) == nt):
        raise ValueError('Error in BTCS:\
                         Argument nt to BTCS should be an integer')
    if not(isinstance(float(d),float) and float(d) > 0):
        raise TypeError('Error in BTCS:\
                        Argument d to BTCS should be a positive float')
    if not(isinstance(phi,np.ndarray)):
        raise TypeError('Error in BTCS:\
                        Argument phi to BTCS should be an array')
    
    nx = len(phi)
     
    # array representing BTCS
    M = np.zeros([nx,nx])
    
    # Zero gradient boundary conditions
    M[0,0] = 1.
    M[0,1] = -1.
    M[-1,-1] = 1.
    M[-1,-2] = -1.
    for i in range(1,nx-1):
        M[i,i-1] = -d
        M[i,i] = 1.+2*d
        M[i,i+1] = -d
    
    # BTCS for all time steps
    for it in range(int(nt)):
        # RHS for zero gradient boundary conditions
        phi[0] = 0
        phi[-1] = 0
        
        phi = la.solve(M,phi)
    
    return phi














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
    there.
    """
    nx = len(phiOld)
    
    # new time-step array for phi
    phi = phiOld.copy()
    
    # FTCS for all time steps
    for it in range(int(nt)):
        
        # initial endpoint for spatial coordinates
        phi[0] = phiOld[0] + d*(phiOld[j] - 2*phiOld[j+1] + phiOld[j+2])
        
        # spatial points, endpoints excluded
        for j in range(1,nx-1):
            phi[j] = phiOld[j] + d*(phiOld[j+1] - 2*phiOld[j] + phiOld[j-1])
       
        # final endpoint for spatial coordinates
        phi[nx-1] = phiOld[0] + d*(phiOld[j-2] - 2*phiOld[j-1] + phiOld[j])
        
        # output to phiOld for the next time-step
        phiOld = phi
        
    return phi
    
def BTCS(phi, d, nt):
    
    """
    Diffusion of profile in phi using BTCS using non-dimensional diffusion
    coefficient d, assuming fixed value boundary conditions.
    """
    nx = len(phi)
     
    # array representing BTCS
    M = np.zeros([nx,nx])
    
    # Zero gradient boundary conditions
    M[0,0] = 1
    M[0,1] = -1
    M[-1,-1] = 1
    M[-1,-2] = -1
    for i in range(1,nx-1):
        M[i,i-1] = -d
        M[i,i] = 1+2*d
        M[i,i+1] = -d
    
    # BTCS for all time steps
    for it in range(int(nt)):
        # RHS for zero gradient boundary conditions
        phi[0] = 0
        phi[-1] = 0
        
        phi = la.solve(M,phi)
    
    return phi














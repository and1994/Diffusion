#!/usr/bin/python

# Outer code for setting up the diffusion problem on a uniform
# grid and calling the function to perform the diffusion and plot.

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
execfile("diffusionSchemes.py")
execfile("diagnostics.py")
execfile("initialConditions.py")

def main():
    """
    Diffuse a squareWave between squareWaveMin and squareWaveMax on a domain
    between x = xmin and x = xmax split over nx spatial steps with diffusion
    coefficient K, time step dt for nt time steps.
    """
    # parameters
    xmin = 0.
    xmax = 1.
    nx = 21
    nt = 20
    dt = 0.1
    K = 2e-3
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    
    # derived parameters
    dx = (xmax - xmin)/(nx-1)
    d = K*dt/dx**2   # non-dimensional diffusion coefficient
    print("non-dimensional diffusion coefficient = ", d)
    print("dx = ", dx, " dt = ", dt, " nt = ", nt)
    print("end time = ", nt*dt)
    
    # spatial points for plotting and for defining initial conditions
    x = np.zeros(nx)
    for j in range(nx):
        x[j] = xmin + j*dx
    print('x = ', x)
    
    # initial conditions
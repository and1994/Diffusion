#!/usr/bin/python

# Outer code for setting up the diffusion problem on a uniform
# grid and calling the function to perform the diffusion and plot.

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
runfile("diffusionSchemes.py")
runfile("diagnostics.py")
runfile("initialConditions.py")

def main(xmin = 0., xmax = 1., nx = 41, nt = 40, dt = 0.1, K = 1e-3, \
         squareWaveMin = 0.4, squareWaveMax = 0.6, name_fig='figure'):
    """
    Diffuse a squareWave between squareWaveMin and squareWaveMax on a domain
    between x = xmin and x = xmax split over nx spatial steps with diffusion
    coefficient K, time step dt for nt time steps.
    """
    # default parameters set in the function arguments
    
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
    phiOld = squareWave(x, squareWaveMin, squareWaveMax)
    
    # analytic solution (of square wave profile in an infinite domain)
    phiAnalytic = analyticErf(x, K*dt*nt, squareWaveMin, squareWaveMax)
    
    # diffusion using FTCS and BTCS
    phiFTCS = FTCS(phiOld.copy(), d, nt)
    phiBTCS = BTCS(phiOld.copy(), d, nt)
    
    # calculate and print out error norms
    print("FTCS L2 error norm = ", L2ErrorNorm(phiFTCS, phiAnalytic))
    print("BTCS L2 error norm = ", L2ErrorNorm(phiBTCS, phiAnalytic))
    
    # plot the solutions
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', linestyle='--', \
             linewidth=2)
    plt.plot(x, phiFTCS, label='FTCS', color='blue')
    plt.plot(x, phiBTCS, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.xlabel('$x$')
    plt.savefig('Plots/'+name_fig+'.pdf')
    
    # plot the errors
    plt.figure(2)
    plt.clf()
    plt.ion()
    plt.plot(x, phiFTCS - phiAnalytic, label='Error FTCS', color='blue')
    plt.plot(x, phiBTCS - phiAnalytic, label='Error BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.05,0.05])
    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.xlabel('$x$')
    plt.savefig('Plots/'+name_fig+'_errors.pdf')
    
main(name_fig='attempt')
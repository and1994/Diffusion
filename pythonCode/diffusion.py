#!/usr/bin/python

# Outer code for setting up the diffusion problem on a uniform
# grid and calling the function to perform the diffusion and plot.

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application (substitute with execfile if supported)
runfile("diffusionSchemes.py")
runfile("diagnostics.py")
runfile("initialConditions.py")

def main(xmin = 0., xmax = 1., nx = 41, d = 0.16, K = 1e-3, \
         squareWaveMin = 0.4, squareWaveMax = 0.6, name_fig='Question3'):
    """
    Diffuse a squareWave between squareWaveMin and squareWaveMax on a domain
    between x = xmin and x = xmax split over nx spatial steps with diffusion
    coefficient K, time step dt for nt time steps.
    """
    # default parameters set in the function arguments
    
    # derived parameters
    dx = (xmax - xmin)/(nx-1)
    dt = d*dx**2/K   # time step imposing the desired value of d
    nt = int(4.0/dt) # imposing same duration for every run
    #print("non-dimensional diffusion coefficient = ", d)
    #print("dx = ", dx, " dt = ", dt, " nt = ", nt)
    #print("end time = ", nt*dt)
    
    # spatial points for plotting and for defining initial conditions
    x = np.zeros(nx)
    for j in range(nx):
        x[j] = xmin + j*dx
    #print('x = ', x)
    
    # initial conditions
    phiOld = squareWave(x, squareWaveMin, squareWaveMax)
    
    # analytic solution (of square wave profile in an infinite domain)
    phiAnalytic = analyticErf(x, K*dt*nt, squareWaveMin, squareWaveMax)
    
    # diffusion using FTCS and BTCS
    phiFTCS = FTCS(phiOld.copy(), d, nt)
    phiBTCS = BTCS(phiOld.copy(), d, nt)
    
    # calculate and print out error norms
    L2errFTCS = L2ErrorNorm(phiFTCS, phiAnalytic)
    L2errBTCS = L2ErrorNorm(phiBTCS, phiAnalytic)
    #print("FTCS L2 error norm = ", L2errFTCS)
    #print("BTCS L2 error norm = ", L2errBTCS)
    
    """
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
    plt.title("t = {:.1f}, d = {:.2f}".format(nt*dt, d))
    plt.savefig('Plots/' + name_fig + '(t=' + str(int(nt*dt)) + ').pdf')
    
    # plot the errors
    plt.figure(2)
    plt.clf()
    plt.ion()
    
    # defining the error vectors that are used in the graph and for evaluating
    # the extremes (m) on the y-axis for plotting
    errorFTCS = phiFTCS - phiAnalytic
    errorBTCS = phiBTCS - phiAnalytic
    m = max(abs(errorFTCS).max(), abs(errorBTCS).max())
    
    plt.plot(x, errorFTCS, label='Error FTCS', color='blue')
    plt.plot(x, errorBTCS, label='Error BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-m,m])
    plt.legend()
    plt.xlabel('$x$')
    plt.title("t = {:.1f}, d = {:.2f}".format(nt*dt, d))
    plt.savefig('Plots/' + name_fig + '(t=' + str(int(nt*dt)) + ')_errors.pdf')
    """
    
    return dx, L2errFTCS, L2errBTCS

def nrms_error_graph(N):
    vector = np.zeros((N,3))
    for it in range(N):
        vector[it,:] = main(nx=21+it*50,d=0.16)
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(3)
    plt.clf()
    plt.ion()
    plt.loglog(vector[:,0].transpose(), vector[:,1].transpose(),'bx',\
               label='FTCS')
    plt.loglog(vector[:,0].transpose(), vector[:,2].transpose(),'ro',\
               label='BTCS')
    plt.loglog(vector[:,0].transpose(),8*(vector[:,0].transpose())**2)
    plt.legend()
    plt.xlabel('$\Delta x$')
    plt.ylabel('$\ell_2$ Norm Error')
    plt.title('Log-log plot of $\ell_2$ Norm Error')
    plt.savefig('Plots/L2errorPlot{}.pdf'.format(int(N)))
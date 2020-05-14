'''
Author:  Kellie McGuire     kellie@kelliejensen.com

Solves the equation of motion for a damped driven pendulum using scipy's
ODEint solver and generates bifurcation diagrams and Poincare sections
'''

import sys
import numpy as np
import scipy.optimize
from scipy.integrate import odeint  ##For solving ODEs
from DDP_plots import bifur_plt, poincare_plt  ##File for generating plots



def main(phi0, omega0, type):

    ###The equation of motion is a second-order ODE
    ###To solve, we write as two first-order ODEs and solve the system
    def f(y0, t, params):
        phi, omega = y0  #Unpack initial conditions
        Q, gamma, Omega, w_0 = params  #Unpack parameters
        derivs = [ omega,
                  -omega/Q  - w_0*w_0*np.sin(phi) + gamma*w_0*w_0*np.cos(Omega*t) ]
        return derivs


    ## DDP parameters
    gamma = 1.5  #Drive strength
    Omega = 2.*np.pi  #Drive frequency; choose 2.*np.pi to make period = 1
    w_0 = 1.5*Omega  #Natural frequency of pendulum
    Q = 2./(w_0)  #Inverse damping
    params = [Q, gamma, Omega, w_0]  #Bundle parameters for solver

    y0 = [phi0, omega0]   #Bundle initial conditions for solver


#############################################################
#      Function for generating data for bifurcation diagrams
#############################################################
    def bifurcation(phi0, omega0):
        ## Time array for solution
        tStop = 600.
        tInc = 0.001
        t = np.arange(0., tStop, tInc)


        # Specify range of drive periods to plot
        period = 2.*np.pi/Omega
        bifur_start = int(500./tInc) #Select bifurcation data afer transients have died out
        bifur_stop = int(600./tInc)
        bifur_Inc = int(period/tInc) #Select data points that are one drive period apart

        gamma_range = np.arange(1.06,1.3,.0001)  #Select range of drive strengths


        # Initialize array to hold data
        bifur_slice = np.empty((len(gamma_range),100,3)) #2nd dim. = number of drive cycles collecting over

        # Calculate phi(t) and omega(t) for range of gamma values
        for i in range(0, len(gamma_range)):
            params[1] = gamma_range[i]
            psoln = odeint(f, y0, t, args=(params,))
            #Select values one drive period apart
            bifur_slice[i,:,0] = psoln[bifur_start:bifur_stop:bifur_Inc, 0]
            bifur_slice[i,:,1] = psoln[bifur_start:bifur_stop:bifur_Inc, 1]
            bifur_slice[i,:,2] = gamma_range[i]
            print(i)

        # Write bifurcation data to files
        np.savetxt(f"Bifur_{phi0}_{omega0}_phis.txt", bifur_slice[:,:,0])
        np.savetxt(f"Bifur_{phi0}_{omega0}_omegas.txt", bifur_slice[:,:,1])
        np.savetxt(f"Bifur_{phi0}_{omega0}_gammas.txt", bifur_slice[:,:,2])

        # Plot phi and omega as a function of gamma (drive strength)
        bifur_plt(phi0, omega0)


#############################################################
#      Function for generating data for Poincare sections
#############################################################
    def poincare(phi0, omega0, gamma):
        ## Time array for solution
        tStop = 60000.
        tInc = 0.001
        t = np.arange(0., tStop, tInc)

        #Specify range of drive periods to plot
        period = 2.*np.pi/Omega
        Poincare_start = int(500./tInc)  #Select data afer transients have died out
        Poincare_stop = int(tStop/tInc)
        Poincare_Inc = int(period/tInc)  #Select data points that are one drive period apart

        Poincare_slice = np.empty((int(tStop)-500,2)) #1st dim. = number of drive cycles collecting over

        #Call solver function
        psoln = odeint(f, y0, t, args=(params,))

        #Select data points one drive period apart
        Poincare_slice[:,0] = psoln[Poincare_start:Poincare_stop:Poincare_Inc, 0]
        Poincare_slice[:,1] = psoln[Poincare_start:Poincare_stop:Poincare_Inc, 1]

        #Save data to file
        np.savetxt(f"Poincare_{phi0}_{omega0}_{gamma}.txt", Poincare_slice)

        #Plot Poincare section for
        poincare_plt(phi0, omega0, gamma)

    if type=="P":
        poincare(phi0, omega0, gamma)
    elif type=="B":
        bifurcation(phi0, omega0)
    else:
        exit(1)



if __name__ == "__main__":
    if len(sys.argv) < 3:
        phi0 = input("Enter initial angular displacement (in radians)  ")
        omega0 = input("Enter intial angular velocity (in radians per second)")
        type = input("Enter 'B' for bifurcation diagram or 'P' for Poincare section")
    else:
        phi0 = sys.argv[1]
        omega0 = sys.argv[2]
        type = sys.argv[3]
    main(phi0, omega0, type)

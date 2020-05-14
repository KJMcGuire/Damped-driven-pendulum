'''
Author:  Kellie McGuire     kellie@kelliejensen.com

Generates bifurcation diagrams and Poincare sections for the damped driven
oscillator using numerical solution generated from odeint (DDP_v3.py)
'''

import matplotlib
matplotlib.use('pdf')  ##Use non-GUI backend
import matplotlib.pyplot as plt
import numpy as np
import math

def bifur_plt(phi0, omega0):
    phis = np.genfromtxt(f"Bifur_{phi0}_{omega0}_phis.txt")
    omegas = np.genfromtxt(f"Bifur_{phi0}_{omega0}_omegas.txt")
    gammas = np.genfromtxt(f"Bifur_{phi0}_{omega0}_gammas.txt")

    #Enforce pi >= phi >= -pi
    for i in range(0, len(phis[:,0])):
        for j in range (0, len(phis[0,:])):
            if phis[i,j] > np.pi:
                phis[i,j] = phis[i,j] - 2*np.pi*math.floor(phis[i,j]/(2*np.pi))
                if phis[i,j] > np.pi:
                    phis[i,j] = phis[i,j]-2*np.pi
            elif phis[i,j] < -np.pi:
                phis[i,j] = phis[i,j] + 2*np.pi*math.floor(abs(phis[i,j])/(2*np.pi))
                if phis[i,j] < -np.pi:
                    phis[i,j] = phis[i,j] + 2*np.pi

    fig1, axs1 = plt.subplots(1, figsize=(8,6))
    fig1.suptitle("Bifurcation Diagram for Damped Driven Pendulum \n $\\phi(0)$ = {}, $\\dot{{\\phi}}(0)$ = {}".format(phi0,omega0))
    axs1.scatter(gammas, phis, s=0.5, c='#000000', marker='.')
    axs1.text(0.23, 0.95,'t = 500, 501, 502, ... 600',
      horizontalalignment='center',
      verticalalignment='top',
      transform = axs1.transAxes)
    axs1.set(xlabel= "drive strength", ylabel = "$\\phi(t)$")
    #axs1.set_xlim([1.05, 1.1])
    #axs1.set_ylim([-1.57, 1.57])
    fig1.savefig(f"DDP_{phi0}_{omega0}_phis.png")

    fig2, axs2 = plt.subplots(1, figsize=(8,6))
    fig2.suptitle("Bifurcation Diagram for Damped Driven Pendulum \n $\\phi(0)$ = {}, $\\dot{{\\phi}}(0)$ = {}".format(phi0,omega0))
    axs2.scatter(gammas, omegas, s=0.5, c='#000000', marker='.')
    axs2.text(0.23, 0.95,'t = 500, 501, 502, ... 600',
      horizontalalignment='center',
      verticalalignment='top',
      transform = axs2.transAxes)
    axs2.set(xlabel= "drive strength", ylabel = "$\\dot{\\phi}(t)$")
    fig2.savefig(f"DDP_{phi0}_{omega0}_omegas.png")


def poincare_plt(phi0, omega0, gamma):
    data = np.genfromtxt(f"Poincare_{phi0}_{omega0}_{gamma}.txt")
    phis = data[:,0]
    omegas = data[:,1]

    #Enforce pi >= phi >= -pi
    for i in range(0,len(phis)):
        if phis[i] > np.pi:
            phis[i] = phis[i] - 2*np.pi*math.floor(phis[i]/(2*np.pi))
            if phis[i] > np.pi:
                phis[i] = phis[i]-2*np.pi
        elif phis[i] < -np.pi:
            phis[i] = phis[i] + 2*np.pi*math.floor(abs(phis[i])/(2*np.pi))
            if phis[i] < -np.pi:
                phis[i] = phis[i] + 2*np.pi

    fig1, axs1 = plt.subplots(1, figsize=(8,6))
    fig1.suptitle("Poincar$\\acute{{e}}$ Section for Damped Driven Pendulum \n $\\phi(0)$ = {}, $\\dot{{\\phi}}(0)$ = {}".format(phi0,omega0))
    axs1.scatter(phis, omegas, s=0.1, c='#000000', marker='.')
    axs1.text(0.23, 0.95,'$\\gamma = {}$'.format(gamma),
      horizontalalignment='center',
      verticalalignment='top',
      transform = axs1.transAxes)
    axs1.set(xlabel= "$\\phi(t)$", ylabel = "$\\dot{\\phi}(t)$")
    fig1.savefig(f"Poincare_{phi0}_{omega0}.png")

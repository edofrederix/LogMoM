#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import figStyle as fs
import os, sys

# Prepare figure

fs.prep(plt)

# Properties

data = np.loadtxt('properties.txt')

M0_0 = data[0]
M2_0 = data[1]
M3_0 = data[2]
case = data[3]

M2_0_tilde = M2_0/(M0_0*M3_0**2)**(1.0/3.0)

if case == 0:

    B1 = 1.0

    TAU = 1.0/B1

    def solution(t_tilde):

        M0_tilde = np.exp(t_tilde)
        M2_tilde = M2_0_tilde*np.exp(t_tilde/5.0)

        return M0_tilde, M2_tilde, np.square(M0_tilde)/np.power(M2_tilde,9.0)

elif case == 1:

    B1 = (M0_0/M3_0)**(1.0/3.0)

    TAU = 1.0/B1*(M0_0/M3_0)**(1.0/3.0)

    def solution(t_tilde):

        M0_tilde = np.power(1.0/15.0*np.square(t_tilde) + 2.0/3.0*M2_0_tilde*t_tilde + 1.0, 3.0/2.0)
        M2_tilde = t_tilde/5.0 + M2_0_tilde

        return M0_tilde, M2_tilde, np.square(M0_tilde)/np.power(M2_tilde,9.0)

elif case == 2:

    B1 = M0_0/M3_0

    TAU = M0_0/(B1*M3_0)

    def solution(t_tilde):

        M0_tilde = 1.0+t_tilde
        M2_tilde = np.power(6.0/5.0*(0.5*np.square(t_tilde) + t_tilde) + M2_0_tilde**6, 1.0/6.0)

        return M0_tilde, M2_tilde, np.square(M0_tilde)/np.power(M2_tilde,9.0)

elif case == 3:

    B1 = 1.0
    B2 = (M0_0/M3_0)**(1.0/3.0)

    TAU = 1.0/B1

    ZETA = B2/B1*(M3_0/M0_0)**(1.0/3.0)

    def solution(t_tilde):

        F0 = \
            3.0*((ZETA+M2_0_tilde)/(3.0/5.0-2.0) + ZETA/2.0) \
          - 3.0*np.exp(-2.0/3.0)*(np.exp(1.0/5.0)*(ZETA+M2_0_tilde)/(3.0/5.0-2.0) + ZETA/2.0)

        C = 1.0 - 2.0/3.0*ZETA*F0

        F = \
            3.0*np.exp(-2.0/3.0*t_tilde)*(np.exp(t_tilde/5.0)*(ZETA+M2_0_tilde)/(3.0/5.0-2.0) + ZETA/2.0) \
          - 3.0*np.exp(-2.0/3.0)*(np.exp(1.0/5.0)*(ZETA+M2_0_tilde)/(3.0/5.0-2.0) + ZETA/2.0)

        M0_tilde = np.exp(t_tilde)*np.power(C + 2.0/3.0*ZETA*F, 3.0/2.0)
        M2_tilde = (ZETA + M2_0_tilde)*np.exp(t_tilde/5.0) - ZETA

        return M0_tilde, M2_tilde, np.square(M0_tilde)/np.power(M2_tilde,9.0)

elif case == 4:

    B1 = 1.0
    B2 = M0_0/M3_0

    TAU = 1.0/B1

    KAPPA = B2/B1*M3_0/M0_0

    def solution(t_tilde):

        G0 = \
            ((KAPPA+1.0)/(1.0-6.0/5.0) + KAPPA/(6.0/5.0)) \
          - np.exp(-6.0/5.0)*((KAPPA+1.0)*np.exp(1.0)/(1.0-6.0/5.0) + KAPPA/(6.0/5.0))

        C = M2_0_tilde**6 - 6.0/5.0*KAPPA*G0

        G = \
            np.exp(-6.0/5.0*t_tilde)*((KAPPA+1.0)*np.exp(t_tilde)/(1.0-6.0/5.0) + KAPPA/(6.0/5.0)) \
          - np.exp(-6.0/5.0)*((KAPPA+1.0)*np.exp(1.0)/(1.0-6.0/5.0) + KAPPA/(6.0/5.0))

        M0_tilde = (KAPPA+1.0)*np.exp(t_tilde) - KAPPA
        M2_tilde = np.exp(t_tilde/5.0)*np.power(C + 6.0/5.0*KAPPA*G, 1.0/6.0)

        return M0_tilde, M2_tilde, np.square(M0_tilde)/np.power(M2_tilde,9.0)

else:

    sys.exit("Invalid case")


# Plot

N = np.loadtxt('postProcessing/probes/0/N.bubbles')
A = np.loadtxt('postProcessing/probes/0/A.bubbles')
alpha = np.loadtxt('postProcessing/probes/0/alpha.bubbles')

M0_tilde = N[:,1]*1e6/M0_0
M2_tilde = A[:,1]/np.pi/(M0_0*M3_0**2)**(1.0/3.0)
t_tilde = N[:,0]/TAU

fig = plt.figure('M0')

plt.plot(t_tilde, M0_tilde, label='LogMoM')

fig = plt.figure('M2')

plt.plot(t_tilde, M2_tilde, label='LogMoM')

# Analytical solution

(M0a_tilde, M2a_tilde, M6a_tilde) = solution(t_tilde)

fig = plt.figure('M0')

plt.plot(t_tilde, M0a_tilde, '--', label='analytical')

fig = plt.figure('M2')

plt.plot(t_tilde, M2a_tilde, '--', label='analytical')

# Style/save

fig = plt.figure('M0')

plt.xlabel(r'$\tilde{t}$')
plt.ylabel(r'$\tilde{M_0}$')

fs.post(fig, plt.legend())

plt.savefig('M0.pdf')

fig = plt.figure('M2')

plt.xlabel(r'$\tilde{t}$')
plt.ylabel(r'$\tilde{M_2}$')

fs.post(fig, plt.legend())

plt.savefig('M2.pdf')

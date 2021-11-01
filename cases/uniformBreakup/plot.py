#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import figStyle as fs
import os, sys

figSize = 'small'

# Prepare figure

fs.prep(plt, figSize)

# Properties

P2 = 0.55

data = np.loadtxt('properties.txt')

N0 = data[0]
S0 = data[1]
A0 = data[2]
case = data[3]

S0_tilde = S0/(N0*A0**2)**(1.0/3.0)

Pi2 = 2.0*P2 - 1.0

if case == 0:

    B1 = 1.0

    TAU = 1.0/B1

    def solution(t_tilde):

        N_tilde = np.exp(t_tilde)
        S_tilde = S0_tilde*np.exp(Pi2*t_tilde)

        return N_tilde, S_tilde, np.square(N_tilde)/np.power(S_tilde,9.0)

elif case == 1:

    B1 = (N0/A0)**(1.0/3.0)

    TAU = 1.0/B1*(N0/A0)**(1.0/3.0)

    def solution(t_tilde):

        N_tilde = np.power(1.0/3.0*Pi2*np.square(t_tilde) + 2.0/3.0*S0_tilde*t_tilde + 1.0, 3.0/2.0)
        S_tilde = Pi2*t_tilde + S0_tilde

        return N_tilde, S_tilde, np.square(N_tilde)/np.power(S_tilde,9.0)

elif case == 2:

    B1 = N0/A0

    TAU = N0/(B1*A0)

    def solution(t_tilde):

        N_tilde = 1.0+t_tilde
        S_tilde = np.power(6.0*Pi2*(0.5*np.square(t_tilde) + t_tilde) + S0_tilde**6, 1.0/6.0)

        return N_tilde, S_tilde, np.square(N_tilde)/np.power(S_tilde,9.0)

elif case == 3:

    B1 = 1.0
    B2 = (N0/A0)**(1.0/3.0)

    TAU = 1.0/B1

    ZETA = B2/B1*(A0/N0)**(1.0/3.0)

    def solution(t_tilde):

        F0 = \
            3.0*((ZETA+S0_tilde)/(3.0*Pi2-2.0) + ZETA/2.0) \
          - 3.0*np.exp(-2.0/3.0)*(np.exp(Pi2)*(ZETA+S0_tilde)/(3.0*Pi2-2.0) + ZETA/2.0)

        C = 1.0 - 2.0/3.0*ZETA*F0

        F = \
            3.0*np.exp(-2.0/3.0*t_tilde)*(np.exp(Pi2*t_tilde)*(ZETA+S0_tilde)/(3.0*Pi2-2.0) + ZETA/2.0) \
          - 3.0*np.exp(-2.0/3.0)*(np.exp(Pi2)*(ZETA+S0_tilde)/(3.0*Pi2-2.0) + ZETA/2.0)

        N_tilde = np.exp(t_tilde)*np.power(C + 2.0/3.0*ZETA*F, 3.0/2.0)
        S_tilde = (ZETA + S0_tilde)*np.exp(Pi2*t_tilde) - ZETA

        return N_tilde, S_tilde, np.square(N_tilde)/np.power(S_tilde,9.0)

elif case == 4:

    B1 = 1.0
    B2 = N0/A0

    TAU = 1.0/B1

    KAPPA = B2/B1*A0/N0

    def solution(t_tilde):

        G0 = \
            ((KAPPA+1.0)/(1.0-6.0*Pi2) + KAPPA/(6.0*Pi2)) \
          - np.exp(-6.0*Pi2)*((KAPPA+1.0)*np.exp(1.0)/(1.0-6.0*Pi2) + KAPPA/(6.0*Pi2))

        C = S0_tilde**6 - 6.0*KAPPA*Pi2*G0

        G = \
            np.exp(-6.0*Pi2*t_tilde)*((KAPPA+1.0)*np.exp(t_tilde)/(1.0-6.0*Pi2) + KAPPA/(6.0*Pi2)) \
          - np.exp(-6.0*Pi2)*((KAPPA+1.0)*np.exp(1.0)/(1.0-6.0*Pi2) + KAPPA/(6.0*Pi2))

        N_tilde = (KAPPA+1.0)*np.exp(t_tilde) - KAPPA
        S_tilde = np.exp(Pi2*t_tilde)*np.power(C + 6.0*KAPPA*Pi2*G, 1.0/6.0)

        return N_tilde, S_tilde, np.square(N_tilde)/np.power(S_tilde,9.0)

else:

    sys.exit("Invalid case")


# Plot

lambd = np.loadtxt('postProcessing/probes/0/lambda.bubbles')
kappa = np.loadtxt('postProcessing/probes/0/kappa.bubbles')
alpha = np.loadtxt('postProcessing/probes/0/alpha.bubbles')
beta = np.loadtxt('postProcessing/probes/0/beta.bubbles')

N_tilde = lambd[:,1]*alpha[:,1]*1e6/N0
S_tilde = kappa[:,1]*alpha[:,1]/np.pi/(N0*A0**2)**(1.0/3.0)
B_tilde = beta[:,1]*alpha[:,1]*N0/A0**2
t_tilde = lambd[:,0]/TAU

fig = plt.figure('N');

plt.plot(t_tilde, N_tilde, label='LogMoM')

fig = plt.figure('S');

plt.plot(t_tilde, S_tilde, label='LogMoM')

fig = plt.figure('B');

plt.plot(t_tilde, B_tilde, label='LogMoM')

# Analytical solution

(Na_tilde, Sa_tilde, Ba_tilde) = solution(t_tilde)

fig = plt.figure('N');

plt.plot(t_tilde, Na_tilde, '--', label='analytical')

fig = plt.figure('S');

plt.plot(t_tilde, Sa_tilde, '--', label='analytical')

fig = plt.figure('B');

plt.plot(t_tilde, Ba_tilde, '--', label='analytical')

# Style/save

fig = plt.figure('N');

plt.xlabel(r'$\tilde{t}$')
plt.ylabel(r'$\tilde{N}$')

plt.yscale('log')

fs.post(fig, figSize, plt.legend())

plt.savefig('N.pdf')

fig = plt.figure('S');

plt.xlabel(r'$\tilde{t}$')
plt.ylabel(r'$\tilde{S}$')

plt.yscale('log')

fs.post(fig, figSize, plt.legend())

plt.savefig('S.pdf')

fig = plt.figure('B');

plt.xlabel(r'$\tilde{t}$')
plt.ylabel(r'$\tilde{B}$')

plt.yscale('log')

fs.post(fig, figSize, plt.legend())

plt.savefig('B.pdf')

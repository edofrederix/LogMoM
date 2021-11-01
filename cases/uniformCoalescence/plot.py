#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import figStyle as fs
import os, sys

figSize = 'small'

# Prepare figure

fs.prep(plt, figSize)

# Properties

data = np.loadtxt('properties.txt')

N0 = data[0]
S0 = data[1]
A0 = data[2]
case = data[3]

B0 = N0**2*A0**8/S0**9
B0_tilde = B0*N0/A0**2

if case == 1:

    K1 = 1.0/N0

    TAU = 1.0/(K1*N0)

    def solution(t_tilde):

        N_tilde = 1.0/(1.0+t_tilde)
        B_tilde = B0_tilde + 2.0*t_tilde

        return N_tilde, np.power(np.square(N_tilde)/B_tilde, 1.0/9.0), B_tilde

elif case == 2:

    K1 = 1.0/A0

    TAU = 1.0/(K1*A0)

    def solution(t_tilde):

        N_tilde = np.exp(-t_tilde)
        B_tilde = B0_tilde*np.exp(2.0*t_tilde)

        return N_tilde, np.power(np.square(N_tilde)/B_tilde, 1.0/9.0), B_tilde

elif case == 4:

    K1 = 1.0/N0
    K2 = 1.0/A0

    TAU = 1.0/(K2*A0)

    KAPPA = K1*N0/(K2*A0)

    def solution(t_tilde):

        N_tilde = 1.0/((KAPPA+1.0)*np.exp(t_tilde) - KAPPA)
        B_tilde = (KAPPA + B0_tilde)*np.exp(2.0*t_tilde) - KAPPA

        return N_tilde, np.power(np.square(N_tilde)/B_tilde, 1.0/9.0), B_tilde

elif case == 5:

    K1 = 1.0/N0

    TAU = 1.0/(K1*N0)

    def solution(t_tilde):

        F0 = B0_tilde
        Q0 = 1.0+F0**(1.0/9.0)

        N_tilde = 1.0/(1.0+t_tilde*Q0)
        B_tilde = B0_tilde + 2.0*t_tilde*Q0

        return N_tilde, np.power(np.square(N_tilde)/B_tilde, 1.0/9.0), B_tilde

elif case == 6:

    K1 = 1.0/N0

    TAU = 1.0/(K1*N0)

    def solution(t_tilde):

        F0 = B0_tilde
        Q0 = 1.0 + F0**(1.0/9.0) + 2.0*F0**(-1.0/6.0)*F0**(1.0/18.0)*(1.0+F0**(2.0/9.0))*F0**(1.0/3.0)

        N_tilde = 1.0/(1.0+t_tilde*Q0*F0**(1.0/3.0))
        B_tilde = B0_tilde + 2.0*t_tilde*Q0

        return N_tilde, np.power(np.square(N_tilde)/B_tilde, 1.0/9.0), B_tilde

elif case == 7:

    K1 = 1.0/(N0**(5.0/6.0)*A0**(1.0/6.0))

    TAU = 1.0/(K1*N0**(5.0/6.0)*A0**(1.0/6.0))

    def solution(t_tilde):

        F0 = B0_tilde
        Q0 = F0**(1.0/12.0)*F0**(1.0/72.0)*(1.0+2.0*F0**(1.0/18.0)+F0**(1.0/3.0))/F0**(1.0/12.0)

        N_tilde = 1.0/(1.0+t_tilde*Q0/F0**(1.0/12.0))
        B_tilde = B0_tilde + 2.0*t_tilde*Q0

        return N_tilde, np.power(np.square(N_tilde)/B_tilde, 1.0/9.0), B_tilde


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

#B_tilde = np.square(N_tilde)/np.power(S_tilde,9.0)

fig = plt.figure('N');

plt.plot(t_tilde, N_tilde, label='LogMoM')

fig = plt.figure('S');

plt.plot(t_tilde, S_tilde, label='LogMoM')

fig = plt.figure('B');

plt.plot(t_tilde, B_tilde, label='LogMoM')

# Analytical solution

if case < 5:

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

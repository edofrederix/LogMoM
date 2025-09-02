#!/bin/python

# Note: the analytical solutions are based on closure using the zeroth, third
# and sixth order moments. The numerical solution is based on the zeroth, second
# and third order moments. Thus, the agreement between numerical and analytical
# solutions will not be perfect.

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

M6_0 = M0_0**2*M3_0**8/M2_0**9
M6_0_tilde = M6_0*M0_0/M3_0**2

if case == 1:

    K1 = 1.0/M0_0

    TAU = 1.0/(K1*M0_0)

    def solution(t_tilde):

        M0_tilde = 1.0/(1.0+t_tilde)
        M6_tilde = M6_0_tilde + 2.0*t_tilde

        return M0_tilde, np.power(np.square(M0_tilde)/M6_tilde, 1.0/9.0), M6_tilde

elif case == 2:

    K1 = 1.0/M3_0

    TAU = 1.0/(K1*M3_0)

    def solution(t_tilde):

        M0_tilde = np.exp(-t_tilde)
        M6_tilde = M6_0_tilde*np.exp(2.0*t_tilde)

        return M0_tilde, np.power(np.square(M0_tilde)/M6_tilde, 1.0/9.0), M6_tilde

elif case == 4:

    K1 = 1.0/M0_0
    K2 = 1.0/M3_0

    TAU = 1.0/(K2*M3_0)

    KAPPA = K1*M0_0/(K2*M3_0)

    def solution(t_tilde):

        M0_tilde = 1.0/((KAPPA+1.0)*np.exp(t_tilde) - KAPPA)
        M6_tilde = (KAPPA + M6_0_tilde)*np.exp(2.0*t_tilde) - KAPPA

        return M0_tilde, np.power(np.square(M0_tilde)/M6_tilde, 1.0/9.0), M6_tilde

elif case == 5:

    K1 = 1.0/M0_0

    TAU = 1.0/(K1*M0_0)

    def solution(t_tilde):

        F0 = M6_0_tilde
        Q0 = 1.0+F0**(1.0/9.0)

        M0_tilde = 1.0/(1.0+t_tilde*Q0)
        M6_tilde = M6_0_tilde + 2.0*t_tilde*Q0

        return M0_tilde, np.power(np.square(M0_tilde)/M6_tilde, 1.0/9.0), M6_tilde

elif case == 6:

    K1 = 1.0/M0_0

    TAU = 1.0/(K1*M0_0)

    def solution(t_tilde):

        F0 = M6_0_tilde
        Q0 = 1.0 + F0**(1.0/9.0) + 2.0*F0**(-1.0/6.0)*F0**(1.0/18.0)*(1.0+F0**(2.0/9.0))*F0**(1.0/3.0)

        M0_tilde = 1.0/(1.0+t_tilde*Q0*F0**(1.0/3.0))
        M6_tilde = M6_0_tilde + 2.0*t_tilde*Q0

        return M0_tilde, np.power(np.square(M0_tilde)/M6_tilde, 1.0/9.0), M6_tilde

elif case == 7:

    K1 = 1.0/(M0_0**(5.0/6.0)*M3_0**(1.0/6.0))

    TAU = 1.0/(K1*M0_0**(5.0/6.0)*M3_0**(1.0/6.0))

    def solution(t_tilde):

        F0 = M6_0_tilde
        Q0 = F0**(1.0/12.0)*F0**(1.0/72.0)*(1.0+2.0*F0**(1.0/18.0)+F0**(1.0/3.0))/F0**(1.0/12.0)

        M0_tilde = 1.0/(1.0+t_tilde*Q0/F0**(1.0/12.0))
        M6_tilde = M6_0_tilde + 2.0*t_tilde*Q0

        return M0_tilde, np.power(np.square(M0_tilde)/M6_tilde, 1.0/9.0), M6_tilde


else:

    sys.exit("Invalid case")


# Plot

lamb = np.loadtxt('postProcessing/probes/0/lambda.bubbles')
kappai = np.loadtxt('postProcessing/probes/0/kappai.bubbles')
alpha = np.loadtxt('postProcessing/probes/0/alpha.bubbles')

M0_tilde = lamb[:,1]*alpha[:,1]*1e6/M0_0
M2_tilde = kappai[:,1]*alpha[:,1]/np.pi/(M0_0*M3_0**2)**(1.0/3.0)
t_tilde = lamb[:,0]/TAU

fig = plt.figure('M0')

plt.plot(t_tilde, M0_tilde, label='LogMoM')

fig = plt.figure('M2')

plt.plot(t_tilde, M2_tilde, label='LogMoM')

# Analytical solution

(M0_tilde, M2_tilde, M6_tilde) = solution(t_tilde)

fig = plt.figure('M0')

plt.plot(t_tilde, M0_tilde, '--', label='analytical')

fig = plt.figure('M2')

plt.plot(t_tilde, M2_tilde, '--', label='analytical')

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

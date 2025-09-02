#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import figStyle as fs
import os, sys

# Properties

H = 2675600
rho_g = 0.59565
T_sat = 373.15
T_liq = 368.15
kappa_c = 0.6759894595

# Prepare figure

fs.prep(plt)

# Properties

data = np.loadtxt('properties.txt')

case = data[0]
dsm_0 = data[1]
sigma_0 = data[2]
alpha_0 = data[3]
Nu_0 = data[4]
d_0 = data[5]

A_0 = 6.0*alpha_0/dsm_0
N_0 = 6.0*alpha_0/np.pi/dsm_0**3*np.exp(3.0*sigma_0**2)

M0_0 = N_0
M2_0 = A_0/np.pi
M3_0 = alpha_0*6.0/np.pi

# Plot LogMoM

lamb = np.loadtxt('postProcessing/probes/0/lambda.steam')
kappai = np.loadtxt('postProcessing/probes/0/kappai.steam')
alpha = np.loadtxt('postProcessing/probes/0/alpha.steam')

t = alpha[:,0]

fig = plt.figure('A')
plt.plot(t, kappai[:,1]*alpha[:,1], label='LogMoM')

fig = plt.figure('alpha')
plt.plot(t, alpha[:,1], label='LogMoM')

fig = plt.figure('dsm')
plt.plot(t, 6.0/kappai[:,1], label='LogMoM')

sigma = np.sqrt(
    np.log(
        (36.0*np.pi*lamb[:,1]*1e6)**(1/3)/kappai[:,1]
    )
)

fig = plt.figure('sigma')
plt.plot(t, sigma, label='LogMoM')

# Analytical solutions

F_c = - 2.0*(T_sat - T_liq)/(H*rho_g)

if case == 1:

    # Constant Nusselt number solution

    C = F_c*kappa_c*Nu_0

    X_0 = M3_0**(4.0/3.0)

    t_end = -M2_0/(2*C*M0_0)

    t = np.linspace(0, t_end*0.95, 1024)

    M2 = 2.0*C*M0_0*t + M2_0
    X = 4.0*C*(C*M0_0*np.square(t) + M2_0*t)*np.cbrt(M0_0) + M3_0**(4.0/3.0)
    M3 = np.power(X, 3.0/4.0)

else:

    t = np.linspace(0, 20, 1024)

    C = F_c*kappa_c*Nu_0/d_0**2

    M2 = M2_0*np.exp(2.0*C*t)
    M3 = M3_0*np.exp(3.0*C*t)

N = np.ones(len(t))*M0_0
A = M2*np.pi
alpha = M3*np.pi/6.0
dsm = M3/M2

lamb = N/np.maximum(alpha, 1e-6)
kappai = A/np.maximum(alpha, 1e-6)

sigma = np.sqrt(
    np.log(
        (36.0*np.pi*lamb)**(1.0/3.0)/kappai
    )
)

fig = plt.figure('A')
plt.plot(t, A, 'k', label='solution')

fig = plt.figure('alpha')
plt.plot(t, alpha, 'k', label='solution')

fig = plt.figure('dsm')
plt.plot(t, dsm, 'k', label='solution')

fig = plt.figure('sigma')
plt.plot(t, sigma, 'k', label='solution')

# Style/save

fig = plt.figure('A')

plt.xlabel(r'$t$')
plt.ylabel(r'$a_i$')

# plt.yscale('log')

fs.post(fig, plt.legend())

plt.savefig('A.pdf')

##

fig = plt.figure('alpha')

plt.xlabel(r'$t$')
plt.ylabel(r'$\alpha$')

# plt.yscale('log')

fs.post(fig, plt.legend())

plt.savefig('alpha.pdf')

##

fig = plt.figure('dsm')

plt.xlabel(r'$t$')
plt.ylabel(r'$d_{sm}$')

# plt.yscale('log')

fs.post(fig, plt.legend())

plt.savefig('dsm.pdf')

##

fig = plt.figure('sigma')

plt.xlabel(r'$t$')
plt.ylabel(r'$\sigma$')

# plt.yscale('log')

fs.post(fig, plt.legend())

plt.savefig('sigma.pdf')

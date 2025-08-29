#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import figStyle as fs
import os, sys

# Prepare figure

fs.prep(plt)

# Properties

data = np.loadtxt('properties.txt')

model = data[0]
dsm = data[1]
sigma = data[2]
alpha = data[3]

d = np.exp(np.linspace(np.log(1e-4), np.log(0.1), 128))

# Plot initial distribution

fig = plt.figure()

dcm = dsm*np.exp(-2.5*sigma**2)
N = 6.0*alpha/np.pi/np.power(dsm,3.0)*np.exp(3.0*np.square(sigma))

dalphadd = \
    np.pi*N*np.square(d)/(6.0*sigma*np.sqrt(np.pi*2.0)) \
  * np.exp(-np.square(np.log(d/dcm))/(2*np.square(sigma)))

plt.plot(d, dalphadd, label='initial')

# Plot final distribution

if model == 1:

    # LogMoM

    A = np.loadtxt('postProcessing/probes/0/A.air')[-1,1]
    alpha = np.loadtxt('postProcessing/probes/0/alpha.air')[-1,1]
    N = np.loadtxt('postProcessing/probes/0/N.air')[-1,1]*1e6

    dsm = 6.0*alpha/A

    sigma = np.sqrt(np.log(N*np.pi/(6.0*alpha)*np.power(dsm,3.0))/3.0)

    dcm = dsm*np.exp(-2.5*sigma**2)

    dalphadd = \
        np.pi*N*np.square(d)/(6.0*sigma*np.sqrt(np.pi*2.0)) \
      * np.exp(-np.square(np.log(d/dcm))/(2*np.square(sigma)))

    plt.plot([], []) # dummy

    plt.plot(d, dalphadd, label='final')

else:

    # FPT

    nSections = int(data[4])

    d = np.loadtxt('d.txt')
    y = np.loadtxt('y.txt')

    dd = y[1:] - y[:-1]

    labels = ['FPT initial', 'FPT final']

    for j,k in enumerate([0,-1]):

        alpha = np.loadtxt('postProcessing/probes/0/alpha.air')[k,1]

        f = np.zeros(nSections)

        for i in range(0,nSections):
            f[i] = np.loadtxt(
                'postProcessing/probes/0/f' +
                str(i) +
                '.air'
            )[k,1]

        dalphadd = alpha*f/dd

        plt.plot(d, dalphadd, 'o-', mew=0.25, mec='w', label=labels[j])

dsm_sim = np.loadtxt('postProcessing/probes/0/d.air')[-1,1]

print('Sauter mean diameter =', dsm_sim)

f_dsm = np.interp(dsm_sim, d, dalphadd)

plt.plot(dsm_sim, f_dsm, 'ks', mew=0.25, mec='w', label='final dsm')

# Style and save

plt.xlabel(r'$d$')
plt.ylabel(r'$\mathrm{d}\alpha/\mathrm{d}d$ [1/m]')

plt.xscale('log')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('plot.pdf')

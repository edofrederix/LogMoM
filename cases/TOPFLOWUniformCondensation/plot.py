#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import figStyle as fs
import os, sys

Nt = 6

# Prepare figure

fs.prep(plt)

# Properties

mode = np.loadtxt('mode.txt')

##

fig = plt.figure('dsm')
dsm = np.loadtxt('postProcessing/probes/0/d.steam')
plt.plot(dsm[:,0], dsm[:,1])

##

fig = plt.figure('alpha')
alpha = np.loadtxt('postProcessing/probes/0/alpha.steam')
plt.plot(alpha[:,0], alpha[:,1], label=r'$\alpha$')

##

fig = plt.figure('dist')

if mode == 1:

    # LogMoM

    d = np.exp(np.linspace(np.log(1e-4), np.log(0.1), 128))

    t = np.loadtxt('postProcessing/probes/0/A.steam')[:,0]
    A = np.loadtxt('postProcessing/probes/0/A.steam')[:,1]
    alpha = np.loadtxt('postProcessing/probes/0/alpha.steam')[:,1]
    N = np.loadtxt('postProcessing/probes/0/N.steam')[:,1]*1e6

    for i in range(0,Nt):

        j = int(round((len(t)-1)*float(i)/(Nt-1)))

        tj = round(t[j]*10)/10
        Aj = A[j]
        alphaj = alpha[j]
        Nj = N[j]

        dsm = 6.0*alphaj/Aj

        sigma = np.sqrt(np.log(Nj*np.pi/(6.0*alphaj)*np.power(dsm,3.0))/3.0)

        dcm = dsm*np.exp(-2.5*sigma**2)

        dalphadd = \
            np.pi*Nj*np.square(d)/(6.0*sigma*np.sqrt(np.pi*2.0)) \
          * np.exp(-np.square(np.log(d/dcm))/(2*np.square(sigma)))

        plt.plot(d, dalphadd, label=r'$t='+str(tj)+'$')

else:

    # FPT

    nSections = int(np.loadtxt('nSections.txt'))

    d = np.loadtxt('d.txt')
    y = np.loadtxt('y.txt')

    dd = y[1:] - y[:-1]

    labels = ['initial fpt', 'final fpt']

    t = np.loadtxt('postProcessing/probes/0/alpha.steam')[:,0]
    alpha = np.loadtxt('postProcessing/probes/0/alpha.steam')[:,1]

    f =  np.zeros([len(t), nSections])

    for i in range(0,nSections):

        f[:,i] = \
            np.loadtxt(
            'postProcessing/probes/0/f' +
            str(i) +
            '.steam'
        )[:,1]

    for i in range(0,Nt):

        j = int(round((len(t)-1)*float(i)/(Nt-1)))

        tj = round(t[j]*10)/10
        alphaj = alpha[j]
        fj = f[j,:]

        dalphadd = alphaj*fj/dd

        plt.plot(d, dalphadd, 'o-', mew=0.25, mec='w', label='$t='+str(tj)+'$')

# Style and save

fig = plt.figure('dsm')

plt.xlabel(r'$t$')
plt.ylabel(r'$d_{sm}$ [m]')

fs.post(fig)

plt.savefig('dsm.pdf')

##

fig = plt.figure('alpha')

plt.xlabel(r'$t$')
plt.ylabel(r'$\alpha$')

fs.post(fig)

plt.savefig('alpha.pdf')

##

fig = plt.figure('dist')

plt.xlabel(r'$d$')
plt.ylabel(r'$\mathrm{d}\alpha/\mathrm{d}d$ [1/m]')

plt.xscale('log')

fs.post(fig, plt.legend(loc='best'))

plt.savefig('dist.pdf')

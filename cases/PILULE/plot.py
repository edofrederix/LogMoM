#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import figStyle as fs
import os, sys
import re

startIndex = 100

# Prepare figure

fs.prep(plt)

# Read

phaseForceNames = ['forceAir', 'forceWater']

for i,phaseForceName in enumerate(phaseForceNames):

    inputFileName = 'postProcessing/'+phaseForceName+'/0/forces.dat'
    outputFileName = 'postProcessing/'+phaseForceName+'/0/forces_parsed.dat'

    # Remove brackets

    outputFile = open(outputFileName, 'w')

    with open(inputFileName, 'r') as inputFile:

        outputFile.write(re.sub(r"[\([{})\]]", "", inputFile.read()))

    outputFile.close()

    # Read parsed file

    inputFileName = outputFileName

    phaseData = np.loadtxt(inputFileName)

    if i == 0:

        # For the first phase create the array

        data = phaseData

    else:

        # For the second phase sum up the forces and moments of both phases

        nTimes = min(phaseData.shape[0], data.shape[0])

        data = data[startIndex:nTimes,:]
        phaseData = phaseData[startIndex:nTimes,:]

        data[:,1:] = data[:,1:] + phaseData[:,1:]

# Plot

fig = plt.figure('forces')

t = data[:,0]
drag = data[:,3]+data[:,6]
lift = data[:,1]+data[:,4]

plt.plot(t, drag, color='C0', label='Drag force')
plt.plot(t, lift, color='C1', label='Lift force')

# Compute running averages

nInterval = int(max(np.round((nTimes-startIndex)/200),10))

tMean = np.array([])
dragMean = np.array([])
liftMean = np.array([])

upper = nInterval

while upper < nTimes-startIndex-1:

    tMean = np.append(tMean, np.mean(t[upper-nInterval:upper]))
    dragMean = np.append(dragMean, np.mean(drag[:upper]))
    liftMean = np.append(liftMean, np.mean(lift[:upper]))

    upper = upper + nInterval

plt.plot(tMean, dragMean, color='w', lw=1.5)
plt.plot(tMean, liftMean, color='w', lw=1.5)

plt.plot(tMean, dragMean, '--', color='C0', lw=1)
plt.plot(tMean, liftMean, '--', color='C1', lw=1)

plt.plot([], [], '-k', label='Instantaneous')
plt.plot([], [], '--k', label='Running average')

print('Mean drag', dragMean[-1])
print('Mean lift', liftMean[-1])

# Style/save

fig = plt.figure('forces')

plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$F$ [N]')

fs.post(fig, plt.legend())

plt.savefig('forces.pdf')

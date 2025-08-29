# Script to check the realized superficial gas velocity at the outlet

import numpy as np

phiAlpha = np.loadtxt('postProcessing/outletPhiAlpha/0/surfaceFieldValue.dat')
area = 4.142621e-04

print('Jg outlet =', phiAlpha[-1,1]/area)

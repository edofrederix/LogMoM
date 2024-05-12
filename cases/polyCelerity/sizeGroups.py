import sys
import numpy as np
from scipy import special

logarithmic = True

dsm = float(sys.argv[1])
sigma = float(sys.argv[2])
nGroups = int(float(sys.argv[3]))
nSectionsPerGroup = int(float(sys.argv[4]))
alphaInlet = float(sys.argv[5])

nSections = nSectionsPerGroup*nGroups

dcm = dsm*np.exp(-2.5*sigma**2)

if (logarithmic):

    w = 6

    logy = np.linspace(
        np.log(dcm) - w*sigma,
        np.log(dcm) + w*sigma,
        nSections+1
    )

    logx = (logy[:-1] + logy[1:])*0.5

    y = np.exp(logy)
    x = np.exp(logx)

else:

    dMax = 1e-1
    dMin = 1e-6

    y = np.linspace(dMin, dMax, nSections+1)
    x = (y[1:] + y[:-1])/2

f = np.zeros(nSections)

def intLogNorm(x,gamma,mu,sigma):
    return \
        0.5*mu**gamma*np.exp(0.5*sigma**2*gamma**2) \
      * special.erf((np.log(x/mu) - gamma*sigma**2)/(2**0.5*sigma))

# The discrete dsm may not be equal to the imposed dsm. Iteratively adjust the
# count mean diameter so that dsm is satisfied.

for corr in range(0,1):

    M3 = (intLogNorm(y[-1],3,dcm,sigma) - intLogNorm(y[0],3,dcm,sigma))

    for i in range(0,nSections):

        f[i] = (intLogNorm(y[i+1],3,dcm,sigma) - intLogNorm(y[i],3,dcm,sigma))/M3

    v = np.pi/6.0*x**3
    x3 = x**3
    x2 = x**2

    dsm2 = sum(f/v*x3)/sum(f/v*x2)

    dcm = dcm*dsm/dsm2

print('FPT dMin =', x[0])
print('FPT dMax =', x[-1])
print('FPT dsm =', sum(f/v*x3)/sum(f/v*x2))

alpha = np.zeros(nGroups)
sumf = np.zeros(nGroups)

for i in range(0,nGroups):

    j = nSectionsPerGroup*i
    k = nSectionsPerGroup*(i+1)

    alpha[i] = np.sum(f[j:k])*alphaInlet
    sumf[i] = np.sum(f[j:k])

    if sumf[i] <= 0:
        f[j:k] = 1/nSectionsPerGroup
    else:
        f[j:k] = f[j:k]/sumf[i]

np.savetxt('alpha.txt', alpha)
np.savetxt('sumf.txt', sumf)
np.savetxt('d.txt', x)
np.savetxt('f.txt', f)

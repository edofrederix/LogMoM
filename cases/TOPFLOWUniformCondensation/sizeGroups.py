import sys
import numpy as np
from scipy import special

logarithmic = True

dsm = float(sys.argv[1])
sigma = float(sys.argv[2])
nSections = int(float(sys.argv[3]))

dcm = dsm*np.exp(-2.5*sigma**2)

dMax = 1e-1
dMin = 1e-4

if (logarithmic):

    logy = np.linspace(
        np.log(dMin),
        np.log(dMax),
        nSections+1
    )

    logx = (logy[:-1] + logy[1:])*0.5

    y = np.exp(logy)
    x = np.exp(logx)

else:

    y = np.linspace(dMin, dMax, nSections+1)
    x = (y[1:] + y[:-1])/2

f = np.zeros(nSections)

def intLogNorm(x,gamma,mu,sigma):
    return \
        0.5*mu**gamma*np.exp(0.5*sigma**2*gamma**2) \
      * special.erf((np.log(x/mu) - gamma*sigma**2)/(2**0.5*sigma))

# The discrete dsm may not be equal to the imposed dsm. Iteratively adjust the
# count mean diameter so that dsm is satisfied.

for corr in range(0,10):

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

np.savetxt('d.txt', x)
np.savetxt('y.txt', y)
np.savetxt('f.txt', f)

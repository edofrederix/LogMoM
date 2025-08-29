#!/bin/python
import sys

pOut = float(sys.argv[1])
L = float(sys.argv[2])
alphaOut = float(sys.argv[3])

rhoWater = 998
M = 28.9
g = 9.98
R = 8.314
T = 303.15

rhoAir = pOut/(R/(M/1e3)*T)
rho = alphaOut*rhoAir + (1-alphaOut)*rhoWater

pIn = pOut + rho*g*L

print(pIn)

VInWater = (1-alphaOut)
VInAir = alphaOut*(pOut/pIn)**(1.0/3.0)

alphaIn = VInAir/(VInAir + VInWater)

print(alphaIn)

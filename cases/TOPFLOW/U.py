import sys

pOut = float(sys.argv[1])
L = float(sys.argv[2])
alphaIn = float(sys.argv[3])
JgOut = float(sys.argv[4])

rhoWater = 998
M = 28.9
g = 9.98
R = 8.314
Rs = R/(M/1e3)
T = 303.15

rhoAirOut = pOut/(Rs*T)
rhoAirIn = (pOut + g*L*(1-alphaIn)*rhoWater)/(Rs*T - g*L*alphaIn)

JgIn = JgOut*rhoAirOut/rhoAirIn

UgIn = JgIn/alphaIn

print(UgIn)

import gmsh
import math
import sys
import threading

nopopup = False
if '-nopopup' in sys.argv:
    nopopup = True
    sys.argv.remove('-nopopup')

if len(sys.argv) < 2:
    raise Exception('Provide mesh number')

MESH = int(sys.argv[1])

# Flow pipe diameter
D1 = 0.038

# Obstacle diameter
D2 = 0.012

# Inlet length
L_in = 0.5

# Outlet length
L_out = 0.25

# Obstacle length
L_obs = 0.034

# Mesh size
dx = 0.5*D1/MESH

# Stretching factor
stretch = 3.0

# Initialize

gmsh.initialize()

gmsh.option.setNumber("General.NumThreads", 16)
# gmsh.option.setNumber("Geometry.AutoCoherence", 0)
gmsh.model.add("PILULE")

# ---- Points ------------------------------------------------------------------

print("Points")

R1 = D1/2
R2 = D2/2

# Outer pipe

c0 = gmsh.model.geo.addPoint(0,0,-R1,dx)

p0 = []
p0.append(gmsh.model.geo.addPoint(R1,0,-R1,dx))
p0.append(gmsh.model.geo.addPoint(0,R1,-R1,dx))
p0.append(gmsh.model.geo.addPoint(-R1,0,-R1,dx))
p0.append(gmsh.model.geo.addPoint(0,-R1,-R1,dx))

# Obstacle

c0_obs = gmsh.model.geo.addPoint(0,-L_obs/2,0,dx)

p0_obs = []
p0_obs.append(gmsh.model.geo.addPoint(0,  -L_obs/2, R2,dx))
p0_obs.append(gmsh.model.geo.addPoint(R2, -L_obs/2, 0, dx))
p0_obs.append(gmsh.model.geo.addPoint(0,  -L_obs/2,-R2,dx))
p0_obs.append(gmsh.model.geo.addPoint(-R2,-L_obs/2, 0, dx))

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Lines -------------------------------------------------------------------

print("Lines")

l0 = []

for i in range(0,4):
    l0.append(gmsh.model.geo.addCircleArc(p0[i], c0, p0[(i+1) % len(p0)]))

l0_obs = []

for i in range(0,4):
    l0_obs.append(
        gmsh.model.geo.addCircleArc(
            p0_obs[i],
            c0_obs,
            p0_obs[(i+1) % len(p0_obs)]
        )
    )

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Surfaces ----------------------------------------------------------------

print("Surfaces")

# Pipe

s0 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop(l0)
])

e01 = []

for i in range(0,4):
    e01.append(
        gmsh.model.geo.extrude(
            [[1, l0[i]]],
            0,
            0,
            D1,
            [int(round(D1/dx))],
            [1],
            False
        )
    )

s01 = []

for i in range(0,4):
    s01.append(e01[i][1][1])

l1 = []

for i in range(0,4):
    l1.append(e01[i][0][1])

s1 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop(l1)
])

# Obstacle

s0_obs = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop(l0_obs)
])

e01_obs = []

for i in range(0,4):
    e01_obs.append(
        gmsh.model.geo.extrude(
            [[1, l0_obs[i]]],
            0,
            L_obs,
            0,
            [int(round(L_obs/dx))],
            [1],
            False
        )
    )

s01_obs = []

for i in range(0,4):
    s01_obs.append(e01_obs[i][1][1])

l1_obs = []

for i in range(0,4):
    l1_obs.append(e01_obs[i][0][1])

s1_obs = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop(l1_obs)
])

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Volumes -----------------------------------------------------------------

def negate(arr):
    return [-x for x in arr]

# Volume around obstacle

v = gmsh.model.geo.addVolume([
    gmsh.model.geo.addSurfaceLoop(
        [s0]
      + [s1]
      + s01
      + [-s0_obs]
      + [-s1_obs]
      + negate(s01_obs)
    )
])

# Inlet volume

L0_inlet = (L_in - D1/2)

Q_inlet = (L0_inlet - stretch*dx)/(L0_inlet - dx)
N_inlet = int(round(math.log(1.0/stretch)/math.log(Q_inlet) + 1))

e0e = gmsh.model.geo.extrude(
    [[2, s0]],
    0,
    0,
    -L0_inlet,
    [N_inlet],
    [1],
    True
)

v_inlet = e0e[1][1]

s_inlet = e0e[0][1]

s0_inlet = []

for i in range(2,len(e0e)):
    s0_inlet.append(e0e[i][1])

# Outlet volume

L0_outlet = (L_out - D1/2)

Q_outlet = (L0_outlet - stretch*dx)/(L0_outlet - dx)
N_outlet = int(round(math.log(1.0/stretch)/math.log(Q_outlet) + 1))

e1e = gmsh.model.geo.extrude(
    [[2, s1]],
    0,
    0,
    (L_out-D1/2),
    [N_outlet],
    [1],
    True
)

v_outlet = e1e[1][1]

s_outlet = e1e[0][1]

s1_outlet = []

for i in range(2,len(e1e)):
    s1_outlet.append(e1e[i][1])

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Physical groups ---------------------------------------------------------

print("Adding physical groups")

def addPhysicalGroup(dim, arr, name):

    group = gmsh.model.addPhysicalGroup(dim, arr)
    gmsh.model.setPhysicalName(dim, group, name)

# Surfaces

addPhysicalGroup(2, [s_inlet], "inlet")
addPhysicalGroup(2, [s_outlet], "outlet")

addPhysicalGroup(2, s01 + s1_outlet + s0_inlet, "wall_pipe")
addPhysicalGroup(2, s01_obs + [s0_obs, s1_obs], "wall_cylinder")

# Volumes

addPhysicalGroup(3, [v, v_inlet, v_outlet], "internal")

# ---- Finalize ----------------------------------------------------------------

print("Finalizing")

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(3)

gmsh.write("PILULE.msh2")

if not nopopup:
    gmsh.fltk.run()

gmsh.finalize()

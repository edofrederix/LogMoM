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

# Chimney radius
R = 1.5675

# Upper plenum radius
Rp = 3.0

# Bypass channel width
b = 0.0285

# Assembly width
a = 0.13914

# Mesh size
dx = a/MESH

# Cell stretch factor in z
S = 2

# Baffle height
hb = 0.2

# Upper plenum mesh height
hp = dx

# Chimney height starting from the baffles to the upper plenum volume
hc = 9.0 - hb - hp

# Standpipe height
hsp = 0.2

# Standpipe radius
spr = 0.1615/2.0

# Standpipe offsets
spox = 0.34426/2.0 + spr
spoy = 0.1305 + 2*spr

# Tolerance
tol = 1e-12

# Heights
z0 = 0
z1 = z0+hb
z2 = z1+hc
z3 = z2+hp
z4 = z3+hsp

# Initialize

gmsh.initialize()

gmsh.option.setNumber("General.NumThreads", 16)
gmsh.option.setNumber("Geometry.AutoCoherence", 0)
gmsh.model.add("BWRX300")

# ---- Some helper functions ---------------------------------------------------

def addPhysicalGroup(dim, arr, name):

    group = gmsh.model.addPhysicalGroup(dim, arr)
    gmsh.model.setPhysicalName(dim, group, name)

def getPointTagsFromLine(tag):
        up, down = gmsh.model.getAdjacencies(1, tag)
        return down

def getLineTagsFromSurface(tag):
        up, down = gmsh.model.getAdjacencies(2, tag)
        return down

def findLineTag(p1, p2):

    for dim, tag in gmsh.model.getEntities(1):

        tags = getPointTagsFromLine(tag)

        if tags[0] == p1 and tags[1] == p2:
            return tag
        elif tags[0] == p2 and tags[1] == p1:
            return -tag

def findPointTag(p):

    points = gmsh.model.getEntitiesInBoundingBox(
        p[0]-tol,
        p[1]-tol,
        p[2]-tol,
        p[0]+tol,
        p[1]+tol,
        p[2]+tol,
        0
    )

    if len(points) > 0:
        return points[0][1]
    else:
        return -1

def getPointCoordinates(tag):
    return gmsh.model.getValue(0, tag, [])

def get2d(arr2d, i, j):
    if i < len(arr2d):
        if j < len(arr2d[i]):
            return arr2d[i][j]

def permute2d(arr2d):

    perm = []
    for j in range(0,25):

        arr = []

        for i in range(0,len(arr2d)):
            if len(arr2d[i]) > j:
                arr.append(arr2d[i][j])

        if len(arr) > 0:
            perm.append(arr)

    return perm

def negate(arr):
    return [-x for x in arr]

def positive(arr):
    return [abs(x) for x in arr]

def extrude(arr, order, h, dx):

    elements = [(order,v) for v in arr if v > 0]

    # Number of layers
    n = max(int(round(h/dx)),1)

    # Extrude
    data = gmsh.model.geo.extrude(
        elements,
        0,
        0,
        h,
        [n],
        [1],
        True
    )

    ret = []

    c = 0
    for i in range(2,len(data)):
        if data[i][0] == order+1:
            ret.append(data[c:i-1])
            c = i-1

    ret.append(data[c:])

    return ret

# ---- Points ------------------------------------------------------------------

print("Adding points")

# Add grid points

p0_grid_x = []

for i in range(0,25):

    x = max(int(i/2)*(a+b) + (i%2)*b - b/2,0)

    arr = []

    for j in range(0,25):

        y = max(int(j/2)*(a+b) + (j%2)*b - b/2,0)

        if x < y:
            y = min(y, (R**2-x**2)**0.5)
        else:
            x = min(x, (R**2-y**2)**0.5)

        if findPointTag([x,y,z0]) == -1:

            arr.append(gmsh.model.geo.addPoint(x,y,z0,dx))
            gmsh.model.geo.synchronize()

    if len(arr) > 0:
        p0_grid_x.append(arr)

# Permuted grid

p0_grid_y = permute2d(p0_grid_x)

# Standpipe centers

n_pipes = [5, 5, 5, 4, 4, 2]

p3_pipes = []
for i,n in enumerate(n_pipes):

    arr = []

    for j in range(0,n):

        px = i*spox
        py = j*spoy + (i%2)*spoy/2

        if (px**2 + py**2)**0.5 + spr < R:

            arr.append([
                gmsh.model.geo.addPoint(px+spr, py,     z3, dx),
                gmsh.model.geo.addPoint(px,     py-spr, z3, dx),
                gmsh.model.geo.addPoint(px-spr, py,     z3, dx),
                gmsh.model.geo.addPoint(px,     py+spr, z3, dx),
                gmsh.model.geo.addPoint(px,     py,     z3, dx)
            ])

    p3_pipes.append(arr)

# Sync

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Lines -------------------------------------------------------------------

print("Adding lines")

n_lines = [
    19, 19, 19, 19, 19, 19,
    17, 17, 17, 17,
    15, 15, 15, 15,
    13, 13,
    9, 9,
    5, 5,
]

# Vertical grid lines

l0_y_grid_x = []
for i in range(0,len(n_lines)):

    py = p0_grid_x[i]

    arr = []

    for j in range(0, n_lines[i]):

        p1 = getPointCoordinates(py[j])
        p2 = getPointCoordinates(py[j+1])

        # Points must be on the same x-coordinate
        if abs(p1[0] - p2[0]) < tol:
            arr.append(gmsh.model.geo.addLine(py[j], py[j+1]))

    l0_y_grid_x.append(arr)

# Horizontal grid lines

l0_x_grid_y = []
for j in range(0,len(n_lines)):

    px = p0_grid_y[j]

    arr = []

    for i in range(0,n_lines[j]):

        p1 = getPointCoordinates(px[i])
        p2 = getPointCoordinates(px[i+1])

        # Points must be on the same y-coordinate
        if abs(p1[1] - p2[1]) < tol:
            arr.append(gmsh.model.geo.addLine(px[i], px[i+1]))

    l0_x_grid_y.append(arr)

# Arcs

baffles = [2, 6]

p0_O = findPointTag([0,0,0])

arcs0_x = []
cursor_x = p0_grid_x[0][-1]
for baffle in baffles:
    for k in range(0,2):
        dest = p0_grid_x[baffle*2+k][-1]
        arcs0_x.append(gmsh.model.geo.addCircleArc(cursor_x, p0_O, dest))
        cursor_x = dest

arcs0_y = []
cursor_y = p0_grid_y[0][-1]
for baffle in baffles:
    for k in range(0,2):
        dest = p0_grid_y[baffle*2+k][-1]
        arcs0_y.append(gmsh.model.geo.addCircleArc(dest, p0_O, cursor_y))
        cursor_y = dest

arcs0 = []

for arc in arcs0_x:
    arcs0.append(arc)

arcs0.append(gmsh.model.geo.addCircleArc(cursor_x, p0_O, cursor_y))

for arc in reversed(arcs0_y):
    arcs0.append(arc)

# Baffle lines

l0_y_baffles = []
for i in baffles:
    for k in range(0,2):
        l0_y_baffles.append(
            gmsh.model.geo.addLine(
                p0_grid_x[i*2+k][-2],
                p0_grid_x[i*2+k][-1]
            )
        )

l0_x_baffles = []
for j in baffles:
    for k in range(0,2):
        l0_x_baffles.append(
            gmsh.model.geo.addLine(
                p0_grid_y[j*2+k][-2],
                p0_grid_y[j*2+k][-1]
            )
        )

# Side lines

l0_sides = [
    gmsh.model.geo.addLine(
        p0_grid_x[0][-2],
        p0_grid_x[0][-1]
    ),
    gmsh.model.geo.addLine(
        p0_grid_y[0][-2],
        p0_grid_y[0][-1]
    )
]

# Standpipe arcs

l3_pipes = []

for py in p3_pipes:

    arr = []

    for pipe in py:

        arr_q = []
        for k in range(0,4):

            p0 = pipe[4]
            p1 = pipe[k]
            p2 = pipe[(k+1)%4]

            c1 = getPointCoordinates(p1)
            c2 = getPointCoordinates(p2)

            in1 = c1[0] > -tol and c1[1] > -tol
            in2 = c2[0] > -tol and c2[1] > -tol

            if in1 and in2:
                arr_q.append(
                    gmsh.model.geo.addCircleArc(
                        p1,
                        p0,
                        p2
                    )
                )

            elif in1 and not in2:
                arr_q.append(
                    gmsh.model.geo.addLine(
                        p1,
                        p0
                    )
                )

            elif in2 and not in1:
                arr_q.append(
                    gmsh.model.geo.addLine(
                        p0,
                        p2
                    )
                )

        arr.append(arr_q)

    l3_pipes.append(arr)

# Lines between standpipes on x0 and y0 edges (kind of hard-coded)

l3_between_pipes_x0 = [
    gmsh.model.geo.addLine(
        p3_pipes[0][0][3],
        p3_pipes[0][1][1]
    ),
    gmsh.model.geo.addLine(
        p3_pipes[0][1][3],
        p3_pipes[0][2][1]
    ),
    gmsh.model.geo.addLine(
        p3_pipes[0][2][3],
        p3_pipes[0][3][1]
    ),
    gmsh.model.geo.addLine(
        p3_pipes[0][3][3],
        p3_pipes[0][4][1]
    )
]

l3_between_pipes_y0 = [
    gmsh.model.geo.addLine(
        p3_pipes[0][0][0],
        p3_pipes[2][0][2]
    ),
    gmsh.model.geo.addLine(
        p3_pipes[2][0][0],
        p3_pipes[4][0][2]
    )
]

l3_x0 = [
    l3_pipes[0][0][1],
    l3_between_pipes_x0[0],
    l3_pipes[0][1][1],
    l3_pipes[0][1][2],
    l3_between_pipes_x0[1],
    l3_pipes[0][2][1],
    l3_pipes[0][2][2],
    l3_between_pipes_x0[2],
    l3_pipes[0][3][1],
    l3_pipes[0][3][2],
    l3_between_pipes_x0[3],
    l3_pipes[0][4][1],
    l3_pipes[0][4][2]
]

l3_y0 = [
    -l3_pipes[0][0][0],
    l3_between_pipes_y0[0],
    -l3_pipes[2][0][1],
    -l3_pipes[2][0][0],
    l3_between_pipes_y0[1],
    -l3_pipes[4][0][1],
    -l3_pipes[4][0][0]
]

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Surfaces ----------------------------------------------------------------

print("Adding surfaces")

n_surfaces = [
    19, 19, 19, 19, 19,
    17, 17, 17, 17,
    15, 15, 15, 15,
    13, 13,
    9, 9,
    5, 5
]

# Grid surfaces

s0_grid_x = []

for i in range(0,19):

    arr = []

    for j in range(0, n_surfaces[i]):

        ll = get2d(l0_y_grid_x,i,  j)
        lr = get2d(l0_y_grid_x,i+1,j)

        lb = get2d(l0_x_grid_y,j,  i)
        lt = get2d(l0_x_grid_y,j+1,i)

        if ll and lr and lb and lt:
            arr.append(gmsh.model.geo.addPlaneSurface([
                gmsh.model.geo.addCurveLoop([
                    lb,
                    lr,
                    -lt,
                    -ll
                ])
            ]))

    if len(arr) > 0:
        s0_grid_x.append(arr)

s0_grid_y = permute2d(s0_grid_x)

# Flatten

s0_grid = []
for sy in s0_grid_x:
    s0_grid += sy

# Set mesh sizes

gmsh.model.geo.synchronize()

for surf in s0_grid:

    lines = getLineTagsFromSurface(surf)

    for line in lines:

        points = getPointTagsFromLine(line)

        p1 = getPointCoordinates(points[0])
        p2 = getPointCoordinates(points[1])

        dist = 0

        for i in range(0,3):
            dist = dist + (p1[i]-p2[i])**2

        dist = dist**0.5

        n = max(int(round(dist/dx)),1)

        gmsh.model.geo.mesh.setTransfiniteCurve(line, n+1)

    gmsh.model.geo.mesh.setTransfiniteSurface(surf)

# Baffle extension surfaces

gmsh.model.geo.synchronize()

s0_y_baffle_extensions = []

for i in range(0,len(baffles)):

    ll = l0_y_baffles[i*2]
    lr = l0_y_baffles[i*2+1]

    pl = getPointTagsFromLine(ll)
    pr = getPointTagsFromLine(lr)

    lb = findLineTag(pl[0], pr[0])
    lt = findLineTag(pl[1], pr[1])

    s0_y_baffle_extensions.append(
        gmsh.model.geo.addPlaneSurface([
            gmsh.model.geo.addCurveLoop([
                lb,
                lr,
                -lt,
                -ll
            ])
        ])
    )

s0_x_baffle_extensions = []

for i in range(0,len(baffles)):

    lb = l0_x_baffles[i*2]
    lt = l0_x_baffles[i*2+1]

    pb = getPointTagsFromLine(lb)
    pt = getPointTagsFromLine(lt)

    ll = findLineTag(pb[0], pt[0])
    lr = findLineTag(pb[1], pt[1])

    s0_x_baffle_extensions.append(
        gmsh.model.geo.addPlaneSurface([
            gmsh.model.geo.addCurveLoop([
                lb,
                lr,
                -lt,
                -ll
            ])
        ])
    )

s0_baffle_extensions = s0_y_baffle_extensions + s0_x_baffle_extensions

# Side surfaces (kind of hard-coded)

s0_sides = [

    gmsh.model.geo.addPlaneSurface([
        gmsh.model.geo.addCurveLoop(
            [arcs0[0]]
          + negate([l0_y_baffles[0]])
          + negate(l0_x_grid_y[19][:4])
          + [l0_sides[0]]
        )
    ]),

    gmsh.model.geo.addPlaneSurface([
        gmsh.model.geo.addCurveLoop(
            [arcs0[2]]
          + negate([l0_y_baffles[2]])
          + negate(reversed(l0_x_grid_y[15][9:12]))
          + l0_y_grid_x[9][15:17]
          + negate(reversed(l0_x_grid_y[17][5:9]))
          + l0_y_grid_x[5][17:19]
          + [l0_y_baffles[1]]
        )
    ]),

    gmsh.model.geo.addPlaneSurface([
        gmsh.model.geo.addCurveLoop(
            [arcs0[4]]
          + negate([l0_x_baffles[3]])
          + negate(reversed(l0_x_grid_y[13][13:15]))
          + l0_y_grid_x[13][13:15]
          + [l0_y_baffles[3]]
        )
    ]),

    gmsh.model.geo.addPlaneSurface([
        gmsh.model.geo.addCurveLoop(
            [arcs0[6]]
          + negate([l0_x_baffles[1]])
          + negate(reversed(l0_x_grid_y[5][17:19]))
          + l0_y_grid_x[17][5:9]
          + negate(reversed(l0_x_grid_y[9][15:17]))
          + l0_y_grid_x[15][9:12]
          + [l0_x_baffles[2]]
        )
    ]),

    gmsh.model.geo.addPlaneSurface([
        gmsh.model.geo.addCurveLoop(
            [arcs0[8]]
          + negate([l0_sides[1]])
          + l0_y_grid_x[19][:4]
          + [l0_x_baffles[0]]
        )
    ]),
]

# Inlets

inlet_high = [
    [0, 3],
    [0, 2],
    [0, 1]
]

inlet_medium = [
    [3,6],
    [2,6],
    [1,5],
    [0,4],
    [0,3],
    [0,2]
]

inlet_low = [
    [6,9],
    [6,9],
    [5,8],
    [4,8],
    [3,7],
    [2,7],
    [0,6],
    [0,4],
    [0,2]
]

s0_inlet_high = []
for i,inlet in enumerate(inlet_high):
    for j in range(inlet[0], inlet[1]):
        s0_inlet_high.append(s0_grid_x[i*2+1][j*2+1])

s0_inlet_medium = []
for i,inlet in enumerate(inlet_medium):
    for j in range(inlet[0], inlet[1]):
        s0_inlet_medium.append(s0_grid_x[i*2+1][j*2+1])

s0_inlet_low = []
for i,inlet in enumerate(inlet_low):
    for j in range(inlet[0], inlet[1]):
        s0_inlet_low.append(s0_grid_x[i*2+1][j*2+1])

s0_inlets = s0_inlet_high + s0_inlet_medium + s0_inlet_low

# Bypasses

s0_y_bypasses = []
for i,sy in enumerate(s0_grid_x[0::2]):
    if i not in baffles:
        for j,surf in enumerate(sy):
            if j/2 not in baffles:
                s0_y_bypasses.append(surf)

s0_x_bypasses = []
for i,sx in enumerate(s0_grid_y[0::2]):
    if i not in baffles:
        for i,surf in enumerate(sx):
            if i/2 not in baffles:
                s0_x_bypasses.append(surf)

s0_bypasses = list(set(s0_y_bypasses) | set(s0_x_bypasses))

# Baffles (all remaining level 0 grid surfaces + baffle extensions)

s0_baffles = \
    list(set(s0_grid) - set(s0_bypasses) - set(s0_inlets)) \
  + s0_baffle_extensions

# Standpipe surfaces, also compute number of cells already

s3_pipes = []
loop3_pipes = []
nc_pipes = []

for i,py in enumerate(l3_pipes):
    for j,pipe in enumerate(py):

        px = i*spox
        py = j*spoy + (i%2)*spoy/2
        r = (px**2 + py**2)**0.5

        t = (Rp**2 - r**2)**0.5 - (Rp**2 - R**2)**0.5
        t0 = Rp - (Rp**2 - R**2)**0.5

        l = hsp + t0 - t

        loop = gmsh.model.geo.addCurveLoop(negate(reversed(pipe)))
        surf = gmsh.model.geo.addPlaneSurface([loop])

        s3_pipes.append(surf)
        loop3_pipes.append(loop)
        nc_pipes.append(int(round(l/dx/S)))

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Extrude first layer -----------------------------------------------------

print("Extruding first layer")

e01 = extrude(
    s0_baffles + s0_inlets + s0_bypasses + s0_sides,
    2,
    z1 - z0,
    dx*S
)

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# Collect level 1 surfaces

s1 = []
for e in e01:
    s1.append(e[0][1])

# Collect baffle surfaces and volumes

s01_baffles = []
v01_baffles = []
for e in e01[:len(s0_baffles)]:

    s01_baffles.append(e[0][1])
    v01_baffles.append(e[1][1])

    for s in e[2:]:
        s01_baffles.append(s[1])

# ---- Extrude second layer ----------------------------------------------------

print("Extruding second layer")

e12 = extrude(
    s1,
    2,
    z2 - z1,
    dx*S
)

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# Collect surfaces

s2 = []
for e in e12:
    s2.append(e[0][1])

# Collect edges

l2_x0 = []
for pi,pj in zip(p0_grid_x[0][:-1], p0_grid_x[0][1:]):

    p1 = getPointCoordinates(pi)
    p2 = getPointCoordinates(pj)

    p1[2] = z2
    p2[2] = z2

    l2_x0.append(
        findLineTag(
            findPointTag(p1),
            findPointTag(p2)
        )
    )

l2_y0 = []
for pi,pj in zip(p0_grid_y[0][:-1], p0_grid_y[0][1:]):

    p1 = getPointCoordinates(pi)
    p2 = getPointCoordinates(pj)

    p1[2] = z2
    p2[2] = z2

    l2_y0.append(
        findLineTag(
            findPointTag(p1),
            findPointTag(p2)
        )
    )

arcs2 = []
for arc in arcs0:

    pi,pj = getPointTagsFromLine(arc)

    p1 = getPointCoordinates(pi)
    p2 = getPointCoordinates(pj)

    p1[2] = z2
    p2[2] = z2

    arcs2.append(
        findLineTag(
            findPointTag(p1),
            findPointTag(p2)
        )
    )

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Extrude third layer ----------------------------------------------------

print("Extruding third layer")

c2_O = [0,0,z2]
c3_O = [0,0,z3]

p2_O = findPointTag(c2_O)
p3_O = findPointTag(c3_O)

# Arcs

arcs3 = []
for arc2 in arcs2:

    sign = -1 if arc2 < 0 else 1

    p21,p22 = getPointTagsFromLine(abs(arc2))

    c31 = getPointCoordinates(p21)
    c32 = getPointCoordinates(p22)

    c31[2] = z3
    c32[2] = z3

    p31 = findPointTag(c31)
    p32 = findPointTag(c32)

    if p31 < 0:
        p31 = gmsh.model.geo.addPoint(c31[0], c31[1], c31[2], dx)

    if p32 < 0:
        p32 = gmsh.model.geo.addPoint(c32[0], c32[1], c32[2], dx)

    gmsh.model.geo.synchronize()

    arc = gmsh.model.geo.addCircleArc(p31, p3_O, p32)
    arcs3.append(sign*arc)

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# Lines

p2_x0 = findPointTag([0,R,z2])
p2_y0 = findPointTag([R,0,z2])

p3_x0 = findPointTag([0,R,z3])
p3_y0 = findPointTag([R,0,z3])

l23_x0 = gmsh.model.geo.addLine(p2_x0, p3_x0)
l23_y0 = gmsh.model.geo.addLine(p2_y0, p3_y0)

l23_arcs = [l23_x0]
for arc2,arc3 in zip(arcs2[:-1],arcs3[:-1]):

    p2 = getPointTagsFromLine(abs(arc2))
    p3 = getPointTagsFromLine(abs(arc3))

    if arc2 > 0:
        p2 = p2[1]
        p3 = p3[1]
    else:
        p2 = p2[0]
        p3 = p3[0]

    line = gmsh.model.geo.addLine(p2,p3)
    l23_arcs.append(line)

l23_arcs.append(l23_y0)

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# Surfaces

s23_wall = []

for i,arc2,arc3 in zip(range(len(arcs2)), arcs2, arcs3):

    surf = gmsh.model.geo.addSurfaceFilling([
            gmsh.model.geo.addCurveLoop([
                -arc2,
                l23_arcs[i],
                arc3,
                -l23_arcs[i+1]
            ])
        ])

    s23_wall.append(surf)


# ---- Make upper plenum -------------------------------------------------------

print("Making upper plenum")

# Lines

l23_O = gmsh.model.geo.addLine(p2_O, p3_O)

l3_x0.append(gmsh.model.geo.addLine(p3_pipes[0][4][3], p3_x0))
l3_y0.append(gmsh.model.geo.addLine(p3_pipes[4][0][0], p3_y0))

# Surfaces

s23_x0 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop(
        negate(reversed(l2_x0))
      + [l23_O]
      + l3_x0
      + [-l23_x0]
    )
])

s23_y0 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop(
        l2_y0
      + [l23_y0]
      + negate(reversed(l3_y0))
      + [-l23_O]
    )
])

s3_plenum = gmsh.model.geo.addPlaneSurface(
    [
        gmsh.model.geo.addCurveLoop(
            negate(reversed(arcs3))
          + negate(reversed(l3_x0))
          + l3_y0
        )
    ]
  + loop3_pipes
)

# Volume

v23 = gmsh.model.geo.addVolume([
    gmsh.model.geo.addSurfaceLoop(
        s23_wall
      + negate(s2)
      + [s3_plenum, s23_x0, s23_y0]
      + s3_pipes
    )
])

# ---- Extrude fourth layer ----------------------------------------------------

print("Extruding fourth layer")

e34 = []
for i,pipe in enumerate(s3_pipes):

    e = extrude(
        [pipe],
        2,
        z4 - z3,
        (z4 - z3)/nc_pipes[i]
    )

    e34.append(e[0])

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# Collect surfaces

s4_outlets = []
s34_pipe_walls = []

for e in e34:
    s4_outlets.append(e[0][1])

    for surf in e[2:]:

        # If any of the points of the surface is not on the x = 0 or y = 0
        # planes then the surface is indeed a wall

        wall_x = False
        wall_y = False

        lines = getLineTagsFromSurface(surf[1])

        for line in lines:

            points = getPointTagsFromLine(line)

            for point in points:

                c = getPointCoordinates(point)

                if c[0] > tol:
                    wall_x = True

                if c[1] > tol:
                    wall_y = True

        if wall_x or wall_y:
            s34_pipe_walls.append(surf[1])

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

# ---- Physical groups ---------------------------------------------------------

print("Adding physical groups")

# Inlets

addPhysicalGroup(2, s0_inlet_high, "inlet_high")
addPhysicalGroup(2, s0_inlet_medium, "inlet_medium")
addPhysicalGroup(2, s0_inlet_low, "inlet_low")

# Bypasses

addPhysicalGroup(2, s0_bypasses, "inlet_bypasses")

# Baffles

addPhysicalGroup(2, s01_baffles, "wall_baffles")

# Sides

addPhysicalGroup(2, s0_sides, "wall_sides")

# Upper plenum

addPhysicalGroup(2, [s3_plenum], "wall_plenum")

# Outlets

addPhysicalGroup(2, s4_outlets, "outlet")

# Pipe walls

addPhysicalGroup(2, s34_pipe_walls, "wall_pipes")

# Volumes

volumes = []
for v in gmsh.model.getEntities(3):
    volumes.append(v[1])

volumes = list(set(volumes) - set(v01_baffles))

addPhysicalGroup(3, volumes, "internal")

# ---- Finalize ----------------------------------------------------------------

print("Finalizing")

gmsh.model.geo.removeAllDuplicates()
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(3)

gmsh.write("BWRX300.msh2")

if not nopopup:
    gmsh.fltk.run()

gmsh.finalize()

FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

convertToMeters 1;

geometry
{
    cylinder
    {
        type   searchableCylinder;
        point1 (0 0 -100);
        point2 (0 0 100);
        radius VARR1;
    }
    cylinder3
    {
        type   searchableCylinder;
        point1 (0 0 -100);
        point2 (0 0 100);
        radius VARR3;
    }
    obstacle
    {
        type   searchableCylinder;
        point1 (0 -100 0);
        point2 (0  100 0);
        radius VARR2;
    }
}

vertices
(
    // 0
    project (VARR2  VARR1MR2  0) (cylinder obstacle)
    (VARR2  VARFR1    0)
    (VARR2  0         0)

    // 3
    (VARRN  VARHR1 0)
    (VARRM  0      0)

    // 5
    (VARR1SINTHETA0 VARR1COSTHETA0 0)
    (VARR1          0              0)

    // 7
    project (0 VARR1 VARR2) (cylinder obstacle)
    (0  VARFR1   VARR2)
    (0  0        VARR2)

    // 10
    project (VARR2BYSQRT2 VARR1MR2 VARR2BYSQRT2) (cylinder obstacle)
    (VARR2BYSQRT2 VARFR1   VARR2BYSQRT2)
    (VARR2BYSQRT2 0        VARR2BYSQRT2)

    // 13
    (0  VARR1  VARZ2)
    (0  VARFR1 VARZ2)
    (0  0      VARZ2)

    // 16
    (VARHR1 VARHR1 VARZ2)
    (VARFR1 0      VARZ2)

    // 18
    (VARR1BYSQRT2  VARR1BYSQRT2 VARZ2)
    (VARR1  0 VARZ2)

    // 20
    (0  VARR1  VARZ3)
    (0  VARFR1 VARZ3)
    (0  0      VARZ3)

    // 23
    (VARHR1 VARHR1 VARZ3)
    (VARFR1 0      VARZ3)

    // 25
    (VARR1BYSQRT2  VARR1BYSQRT2 VARZ3)
    (VARR1  0 VARZ3)

    // New outer layer points

    // 27
    project (VARR2  VARR3MR2  0) (cylinder3 obstacle)
    (VARR3SINTHETA0 VARR3COSTHETA0 0)
    (VARR3          0              0)

    // 30
    project (0 VARR3 VARR2) (cylinder3 obstacle)
    project (VARR2BYSQRT2 VARR3MR2 VARR2BYSQRT2) (cylinder3 obstacle)

    // 32
    (0  VARR3  VARZ2)
    (VARR3BYSQRT2  VARR3BYSQRT2 VARZ2)
    (VARR3  0 VARZ2)

    // 35
    (0  VARR3  VARZ3)
    (VARR3BYSQRT2  VARR3BYSQRT2 VARZ3)
    (VARR3  0 VARZ3)

    // New points for the slit

    // 38
    (0 VARR3 0)
    (0 VARR1 0)

    // 40
    (VARFR2 VARR3MFR2 0)
    (VARFR2 VARR1MFR2 0)

    // 42
    (0 VARR3 VARFR2)
    (0 VARR1 VARFR2)

    // 44
    (VARHR2 VARR3MFR2 VARHR2)
    (VARHR2 VARR1MFR2 VARHR2)
);

blocks
(
    hex (2 4 3 1 12 17 16 11) (VARNC VARNB VARNB) simpleGrading (1 1 1)
    hex (4 6 5 3 17 19 18 16) (VARNR1 VARNB VARNB) simpleGrading (1 1 1)
    hex (9 12 11 8 15 17 16 14) (VARNB VARNB VARNC) simpleGrading (1 1 1)
    hex (8 11 10 7 14 16 18 13) (VARNB VARNR1 VARNC) simpleGrading (1 1 1)
    hex (11 1 0 10 16 3 5 18) (VARNB VARNR1 VARNC) simpleGrading (1 1 1)

    hex (15 17 16 14 22 24 23 21) (VARNB VARNB VARNL) simpleGrading (1 1 VARGL)
    hex (17 19 18 16 24 26 25 23) (VARNR1 VARNB VARNL) simpleGrading (1 1 VARGL)
    hex (16 18 13 14 23 25 20 21) (VARNR1 VARNB VARNL) simpleGrading (1 1 VARGL)

    hex (6 29 28 5 19 34 33 18) (VARNS VARNB VARNB) simpleGrading (1 1 1)
    hex (5 28 27 0 18 33 31 10) (VARNS VARNC VARNB) simpleGrading (1 1 1)
    hex (7 10 31 30 13 18 33 32) (VARNB VARNS VARNC) simpleGrading (1 1 1)

    hex (19 34 33 18 26 37 36 25) (VARNS VARNB VARNL) simpleGrading (1 1 VARGL)
    hex (18 33 32 13 25 36 35 20) (VARNS VARNB VARNL) simpleGrading (1 1 VARGL)

    hex (39 41 40 38 43 45 44 42) (VARNB VARNS VARNB) simpleGrading (1 1 1)
    hex (41 0 27 40 45 10 31 44) (VARNBH VARNS VARNB) simpleGrading (1 1 1)
    hex (43 45 44 42 7 10 31 30) (VARNB VARNS VARNBH) simpleGrading (1 1 1)
);

edges
(
    project 5 0 (cylinder)
    project 6 5 (cylinder)

    project 18 10 (cylinder)

    project 10 0 (obstacle)
    project 11 1 (obstacle)
    project 12 2 (obstacle)

    project 7 10 (obstacle)
    project 8 11 (obstacle)
    project 9 12 (obstacle)

    project 18 13 (cylinder)
    project 19 18 (cylinder)

    project 25 20 (cylinder)
    project 26 25 (cylinder)

    project 29 28 (cylinder3)
    project 34 33 (cylinder3)
    project 36 35 (cylinder3)
    project 37 36 (cylinder3)

    project 28 27 (cylinder3)
    project 33 32 (cylinder3)

    project 39 41 (cylinder)
    project 41 0 (cylinder)

    project 38 40 (cylinder3)
    project 40 27 (cylinder3)

    project 27 31 (obstacle)
    project 31 30 (obstacle)

    project 31 33 (cylinder3)
);

boundary
(
    wall_cylinder
    {
        type wall;
        faces
        (
            (12 9 8 11)
            (11 8 7 10)
            (11 12 2 1)
            (10 11 1 0)
            (39 43 45 41)
            (45 10 0 41)
            (7 10 45 43)
        );
    }
    wall_pipe
    {
        type wall;
        faces
        (
            (28 33 34 29)
            (33 36 37 34)
            (27 31 33 28)
            (30 32 33 31)
            (32 35 36 33)
            (38 42 44 40)
            (42 30 31 44)
            (44 31 27 40)
        );
    }
    inlet
    {
        type patch;
        faces
        (
            (20 21 23 25)
            (21 22 24 23)
            (23 24 26 25)
            (35 20 25 36)
            (36 25 26 37)
        );
    }
);

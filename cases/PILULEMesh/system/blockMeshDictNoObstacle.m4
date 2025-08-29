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
    project (VARR2  VARR3MR2  0) (cylinder3 obstacle)
    (VARR2  VARFR3    0)
    (VARR2  0         0)

    (VARRN  VARHR3 0)
    (VARRM  0      0)

    (VARR3SINTHETA0 VARR3COSTHETA0 0)
    (VARR3          0              0)

    project (0 VARR3 VARR2) (cylinder3 obstacle)
    (0  VARFR3   VARR2)
    (0  0        VARR2)

    project (VARR2BYSQRT2 VARR3MR2 VARR2BYSQRT2) (cylinder3 obstacle)
    (VARR2BYSQRT2 VARFR3   VARR2BYSQRT2)
    (VARR2BYSQRT2 0        VARR2BYSQRT2)

    (0  VARR3  VARZ2_NO)
    (0  VARFR3 VARZ2_NO)
    (0  0      VARZ2_NO)

    (VARHR3 VARHR3 VARZ2_NO)
    (VARFR3 0      VARZ2_NO)

    (VARR3BYSQRT2  VARR3BYSQRT2 VARZ2_NO)
    (VARR3  0 VARZ2_NO)

    (0  VARR3  VARZ3_NO)
    (0  VARFR3 VARZ3_NO)
    (0  0      VARZ3_NO)

    (VARHR3 VARHR3 VARZ3_NO)
    (VARFR3 0      VARZ3_NO)

    (VARR3BYSQRT2  VARR3BYSQRT2 VARZ3_NO)
    (VARR3  0 VARZ3_NO)
);

blocks
(
    hex (15 17 16 14 22 24 23 21) (VARNB VARNB VARNL_NO) simpleGrading (1 1 1)
    hex (17 19 18 16 24 26 25 23) (VARNR1 VARNB VARNL_NO) simpleGrading (1 1 1)
    hex (16 18 13 14 23 25 20 21) (VARNR1 VARNB VARNL_NO) simpleGrading (1 1 1)
);

edges
(
    project 18 13 (cylinder3)
    project 19 18 (cylinder3)

    project 25 20 (cylinder3)
    project 26 25 (cylinder3)
);

boundary
(
    wall_pipe
    {
        type wall;
        faces
        (
            (18 25 26 19)
            (13 20 25 18)
        );
    }
    inlet
    {
        type patch;
        faces
        (
            (13 14 16 18)
            (14 15 17 16)
            (16 17 19 18)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (20 21 23 25)
            (21 22 24 23)
            (23 24 26 25)
        );
    }
);

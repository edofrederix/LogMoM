FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

convertToMeters 1;

vertices
(
    (0 0 0)
    (0.0960139398 0 0.004192059145)
    (0.0960139398 VARL 0.004192059145)
    (0 VARL 0)
    (0.0960139398 0 -0.004192059145)
    (0.0960139398 VARL -0.004192059145)
    (0.0974072016 0 0.004252890268)
    (0.0974072016 0 -0.004252890268)
    (0.0974072016 VARL -0.004252890268)
    (0.0974072016 VARL 0.004252890268)

);

blocks
(
    hex (0 6 7 0 3 9 8 3) (VARNX 1 VARNZ) simpleGrading (1 1 1)
);

edges
(
    arc 8 9 (0.0975 VARL 0)
    arc 6 7 (0.0975 0 0)
);

patches
(
    wedge wedgeFront
    (
        (0 6 9 3)
    )
    wedge wedgeBack
    (
        (0 3 8 7)
    )
    wall walls
    (
        (6 7 8 9)
    )
    patch inlet
    (
        (0 7 6 0)
    )
    patch outlet
    (
        (3 8 9 3)
    )
    empty axis
    (
        (0 3 3 0)
    )
);
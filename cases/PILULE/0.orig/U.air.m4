FoamFile
{
    version     2.0;
    format      binary;
    class       volVectorField;
    object      U.air;
}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 VARUG);

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           $internalField;
    }

    outlet
    {
        type            pressureInletOutletVelocity;
        phi             phi.air;
        value           $internalField;
    }

    wall_cylinder
    {
        type            VARUAIRCYLINDERBC;
        value           uniform (0 0 0);
    }

    wall_pipe
    {
        type            slip;
    }

    symm
    {
        type            symmetry;
    }
}

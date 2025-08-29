FoamFile
{
    format      ascii;
    class       volScalarField;
    object      alpha.air;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARALPHAG;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           $internalField;
    }

    outlet
    {
        type            zeroGradient;
    }

    "wall.*"
    {
        type            zeroGradient;
    }

    symm
    {
        type            symmetry;
    }
}

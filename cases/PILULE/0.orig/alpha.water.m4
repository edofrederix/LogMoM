FoamFile
{
    format      ascii;
    class       volScalarField;
    object      alpha.water;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARALPHAL;

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

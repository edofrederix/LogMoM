FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T.steam;
}

dimensions          [0 0 0 1 0 0 0];

internalField       uniform 560.15;

boundaryField
{
    inlet_low
    {
        type            fixedValue;
        value           uniform VARTLOW;
    }

    inlet_medium
    {
        type            fixedValue;
        value           $internalField;
    }

    inlet_high
    {
        type            fixedValue;
        value           $internalField;
    }

    inlet_bypasses
    {
        $inlet_low
    }

    outlet
    {
        type            zeroGradient;
    }

    "wall.*"
    {
        type            zeroGradient;
    }

    "symm.*"
    {
        type            symmetry;
    }
}

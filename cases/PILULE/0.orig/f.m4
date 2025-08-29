FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      VARFNAME.air;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARFIELD;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           $internalField;
    }

    outlet
    {
        type            inletOutlet;
        phi             phi.air;
        inletValue      $internalField;
        value           $internalField;
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

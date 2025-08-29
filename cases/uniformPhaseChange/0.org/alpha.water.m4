FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      alpha.water;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARALPHAWATER;

boundaryField
{
    empties
    {
        type            empty;
    }
}

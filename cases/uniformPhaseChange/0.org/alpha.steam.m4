FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      alpha.steam;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARALPHASTEAM;

boundaryField
{
    empties
    {
        type            empty;
    }
}

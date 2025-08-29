FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      VARFNAME.steam;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARFIELD;

boundaryField
{
    empties
    {
        type            empty;
    }
}
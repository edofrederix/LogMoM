FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T.steam;
}

dimensions          [0 0 0 1 0 0 0];

internalField       uniform VARTVAP;

boundaryField
{
    empties
    {
        type            empty;
    }
}

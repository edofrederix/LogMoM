FoamFile
{
    format      ascii;
    class       volScalarField;
    object      alpha.water;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARALPHA;

boundaryField
{
    empties
    {
        type            empty;
    }
}

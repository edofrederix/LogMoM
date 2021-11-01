FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      alpha.bubbles;
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

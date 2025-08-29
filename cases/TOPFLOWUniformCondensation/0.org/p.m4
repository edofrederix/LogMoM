FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}

dimensions          [1 -1 -2 0 0 0 0];

internalField       uniform VARPRESSURE;

boundaryField
{
    empties
    {
        type            empty;
    }
}

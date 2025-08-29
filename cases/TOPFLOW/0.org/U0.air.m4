FoamFile
{
    version     2.0;
    format      binary;
    class       volVectorField;
    object      U0.air;
}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 VARUAIRIN 0);

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           $internalField;
    }
    outlet
    {
        type            pressureInletOutletVelocity;
        phi             phi0.air;
        value           $internalField;
    }
    walls
    {
        type            noSlip;
    }
    wedgeFront
    {
        type            wedge;
    }
    wedgeBack
    {
        type            wedge;
    }
    axis
    {
        type            empty;
    }
}

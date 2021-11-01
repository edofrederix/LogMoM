FoamFile
{
    version     2.0;
    format      binary;
    class       volVectorField;
    object      U.air;
}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 VARUAIR 0);

boundaryField
{
    inlet
    {
        type               fixedValue;
        value              $internalField;
    }
    outlet
    {
        type               pressureInletOutletVelocity;
        phi                phi.air;
        value              $internalField;
    }
    walls
    {
        type               fixedValue;
        value              uniform (0 0 0);
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

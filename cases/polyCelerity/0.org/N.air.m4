FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      N.air;
}

dimensions      [0 -3 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type            inletOutletLogNormal;
        phi             VARFLUX0.air;
        sigma           VARSIGMA;
        dsm             VARDSM;
        value           $internalField;
    }
    outlet
    {
        type            inletOutletLogNormal;
        phi             VARFLUX0.air;
        sigma           VARSIGMA;
        dsm             VARDSM;
        value           $internalField;
    }
    walls
    {
        type            zeroGradient;
    }
    axis
    {
        type            empty;
    }
    wedgeFront
    {
        type            wedge;
    }
    wedgeBack
    {
        type            wedge;
    }
}

FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      lambda.steam;
}

dimensions      [0 -3 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type            inletOutletLogNormal;
        phi             phi.steam;
        sigma           VARSIGMA;
        dsm             VARDSM;
        value           $internalField;
    }
    outlet
    {
        type            inletOutletLogNormal;
        phi             phi.steam;
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

sources
{
    phaseChange
    {
        type            interfacialGrowthMoment;
    }
}

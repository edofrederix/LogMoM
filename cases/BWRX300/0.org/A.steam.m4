FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      A.steam;
}

dimensions      [0 -1 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    "inlet.*"
    {
        type            inletOutletLogNormal;
        phi             phi.steam;
        sigma           VARSIGMA;
        dsm             VARDSM;
        value           $internalField;
    }

    outlet
    {
        type            zeroGradient;
    }

    "wall.*"
    {
        type            zeroGradient;
    }

    "symm.*"
    {
        type            symmetry;
    }
}

sources
{
    phaseChange
    {
        type            interfacialGrowthMoment;
    }
}

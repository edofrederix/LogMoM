FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      A.air;
}

dimensions      [0 -1 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type            inletOutletLogNormal;
        phi             VARFLUX2.air;
        sigma           VARSIGMA;
        dsm             VARDSM;
        value           $internalField;
    }

    outlet
    {
        type            inletOutlet;
        phi             VARFLUX2.air;
        inletValue      $internalField;
        value           $internalField;
    }

    "wall.*"
    {
        type            zeroGradient;
    }
}

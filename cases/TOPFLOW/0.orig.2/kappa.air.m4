FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      kappa.air;
}

dimensions      [0 -1 0 0 0 0 0];

internalField   uniform 0.0;

boundaryField
{
    inlet
    {
        type            inletOutletLogNormal;
        phi             phi.air;
        sigma           VARSIGMA;
        dsm             VARDSM;
        alphaName       air;
        value           $internalField;
    }
    outlet
    {
        type            inletOutlet;
        phi             phi.air;
        inletValue      $internalField;
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

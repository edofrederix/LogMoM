FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.VARPHASENAME;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 1e-5;

boundaryField
{
    wedgeFront
    {
        type            wedge;
    }
    wedgeBack
    {
        type            wedge;
    }
    walls
    {
        type            zeroGradient;
    }
    inlet
    {
        type            TOPFLOWAlphaInlet;
        average         VARALPHAIN;
        case            VARCASE;
        phase           air;
    }
    outlet
    {
        type            inletOutlet;
        phi             phi.VARPHASENAME;
        inletValue      $internalField;
        value           $internalField;
    }
    axis
    {
        type            empty;
    }
}

FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      VARFNAME.VARPHASENAME;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARFIELD;

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
        type            fixedValue;
        value           uniform VARINLET;
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

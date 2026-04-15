FoamFile
{
    format      ascii;
    class       volScalarField;
    object      alpha.water;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARALPHAL;

boundaryField
{
    inlet
    {
        type            oscillatingDisk;
        lambda          0.0095;
        kappa           0.04;
        average         VARALPHAL;
        amplitude       VARAMPLITUDE;
    }

    outlet
    {
        type            zeroGradient;
    }

    "wall.*"
    {
        type            zeroGradient;
    }

    symm
    {
        type            symmetry;
    }
}

FoamFile
{
    format      ascii;
    class       volScalarField;
    object      alpha.air;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform VARALPHAG;

boundaryField
{
    inlet
    {
        type            oscillatingDisk;
        lambda          0.0095;
        kappa           0.04;
        average         VARALPHAG;
        amplitude       -VARAMPLITUDE;
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

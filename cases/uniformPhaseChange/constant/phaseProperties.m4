FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      phaseProperties;
}

phases (steam water);

water
{
    type            purePhaseModel;

    diameterModel   constant;

    constantCoeffs
    {
        d           1e-3;
    }

    residualAlpha   1e-12;
}

steam
{
    type            purePhaseModel;

    diameterModel   LogMoM;

    LogMoMCoeffs
    {
        dMin            1e-4;
        dMax            1e-1;
        sigmaMax        1.0;

        continuousPhase water;

        // To obtain the proper analytical solutions, we must disable the
        // realizability mode. This means that the solution may become
        // unrealizable at some point in time. However, since the moments are
        // limited, some solution will be produced anyway, whereas the
        // analytical solution will no longer exist.

        realizableShrinkage false;

        sources
        ();
    }

    residualAlpha   1e-12;
}

blending
{
    default
    {
        type    continuous;
        phase   water;
    }
}

surfaceTension
{
    steam_water
    {
        type constant;
        sigma 0.072;
    }
}

interfaceCompression
{}

aspectRatio
{}

FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      phaseProperties;
}

phases (water bubbles);

water
{
    type            pureIsothermalPhaseModel;

    diameterModel   constant;

    constantCoeffs
    {
        d           1e-3;
    }

    residualAlpha   1e-6;
}

bubbles
{
    type            pureIsothermalPhaseModel;

    diameterModel   LogMoM;

    LogMoMCoeffs
    {
        dMax            100.0;
        dMin            1e-16;

        continuousPhase water;

        sources
        (
            breakup
            {
                type            polynomial;
                B               (VARBREAKRATE1 VARBREAKRATE2);
                p               (VARR1 VARR2);
                GaussHermite    VARNGH;
                GaussLegendre   VARNGH;
            }
        );
    }

    residualAlpha   1e-6;
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
{}

interfaceCompression
{}

aspectRatio
{}

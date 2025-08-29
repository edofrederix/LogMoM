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
            coalescence
            {
                type            polynomial;
                K               (VARCOARATE1 VARCOARATE2 VARCOARATE3 VARCOARATE4);
                p               (VARP1 VARP2 VARP3 VARP4);
                q               (VARQ1 VARQ2 VARQ3 VARQ4);
                GaussHermite    VARNGH;
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

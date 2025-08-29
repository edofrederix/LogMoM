FoamFile
{
    format      ascii;
    class       dictionary;
    object      phaseProperties;
}

phases (air water);

air
{
    type            pureIsothermalPhaseModel;

    diameterModel   LogMoM;

    LogMoMCoeffs
    {
        dMax            1e-1;
        dMin            1e-4;

        continuousPhase water;

        sources
        (
            coalescence
            {
                type            VARCOALESCENCEMODEL;
                GaussHermite    3;

                turbulence      true;
                buoyancy        false;
                laminarShear    false;
            }

            breakup
            {
                type            VARBREAKUPMODEL;
                GaussHermite    3;
                GaussLegendre   5;
            }
        );
    }

    residualAlpha   1e-6;
}

water
{
    type            pureIsothermalPhaseModel;

    diameterModel   constant;

    constantCoeffs
    {
        d           1e-4;
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
{
    air_water
    {
        type    constant;
        sigma   0.07;
    }
}

interfaceCompression
{}

aspectRatio
{}

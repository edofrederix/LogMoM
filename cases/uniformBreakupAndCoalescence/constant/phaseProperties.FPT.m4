FoamFile
{
    format      ascii;
    class       dictionary;
    object      phaseProperties;
}

phases (air water);

populationBalances (bubbles);

air
{
    type            pureIsothermalPhaseModel;
    diameterModel   velocityGroup;

    velocityGroupCoeffs
    {
        populationBalance    bubbles;

        shapeModel           spherical;

        sizeGroups
        (
            #include "FPT/sizeGroups"
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

populationBalanceCoeffs
{
    bubbles
    {
        continuousPhase water;

        coalescenceModels
        (
            VARCOALESCENCEMODEL
            {
                turbulence      true;
                buoyancy        false;
                laminarShear    false;
            }
        );

        binaryBreakupModels
        (
            VARBREAKUPMODEL
            {}
        );

        breakupModels
        ();
    }
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

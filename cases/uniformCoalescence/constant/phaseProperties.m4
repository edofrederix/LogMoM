FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      phaseProperties;
}

type basicMultiphaseSystem;

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

    diameterModel   threeMomentLogNormal;

    threeMomentLogNormalCoeffs
    {
        dMax            100.0;
        dMin            1e-16;

        continuousPhase water;

        p               3;
        q               2;

        closingMoment   squaredVolume;

        coalescence
        {
            active          true;

            efficiencyType  polynomial;
            frequencyType   constant;

            polynomialEfficiencyCoeffs
            {
                K   (VARCOARATE1 VARCOARATE2 VARCOARATE3 VARCOARATE4);
                p   (VARP1 VARP2 VARP3 VARP4);
                q   (VARQ1 VARQ2 VARQ3 VARQ4);
            }

            constantFrequencyCoeffs
            {
                K   1.0;
            }

            GaussHermite    VARNGH;
        }

        breakup
        {
            active  false;
        }
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
();

interfaceCompression
();

aspectRatio
();

drag
(
    (bubbles in water)
    {
        type            SchillerNaumann;
        residualRe      1e-3;
    }
);

virtualMass
();

heatTransfer
();

phaseTransfer
();

lift
();

wallLubrication
();

turbulentDispersion
();

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

        closingMoment   interfacialArea;

        coalescence
        {
            active      false;
        }

        breakup
        {
            active      true;

            type        polynomial;
            B           (VARBREAKRATE1 VARBREAKRATE2);
            r0          (VARR1 VARR2);

            GaussHermite    VARNGH;
            GaussLegendre   VARNGH;
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

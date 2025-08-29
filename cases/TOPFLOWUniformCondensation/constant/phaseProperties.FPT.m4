FoamFile
{
    format      ascii;
    class       dictionary;
    object      phaseProperties;
}

phases (steam water);

populationBalances (bubbles);

steam
{
    type            purePhaseModel;
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
    type            purePhaseModel;

    diameterModel   constant;

    constantCoeffs
    {
        d           1e-3;
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
            LehrMilliesMewes
            {}
        );

        binaryBreakupModels
        (
            LehrMilliesMewes
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
        type            linear;
        minFullyContinuousAlpha.steam 0.7;
        minPartlyContinuousAlpha.steam 0.3;
        minFullyContinuousAlpha.water 0.7;
        minPartlyContinuousAlpha.water 0.3;
    }

    drag
    {
        type            linear;
        minFullyContinuousAlpha.steam 0.7;
        minPartlyContinuousAlpha.steam 0.5;
        minFullyContinuousAlpha.water 0.7;
        minPartlyContinuousAlpha.water 0.5;
    }
}

surfaceTension
{
    steam_water
    {
        type    constant;
        sigma   VARSTEN;
    }
}

interfaceCompression
{}

aspectRatio
{
    steam_dispersedIn_water
    {
        type Wellek;
    }
}

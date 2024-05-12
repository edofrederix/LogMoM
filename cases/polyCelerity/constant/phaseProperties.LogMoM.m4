FoamFile
{
    format      ascii;
    class       dictionary;
    object      phaseProperties;
}

type    VARPHASESYSTEM;

phases (air water);

air
{
    type            VARPHASEMODEL;

    diameterModel   threeMomentLogNormal;

    threeMomentLogNormalCoeffs
    {
        dMax            5e-2;
        dMin            1e-6;

        p               3;
        q               2;

        continuousPhase water;

        coalescence
        {
            active      false;
        }

        breakup
        {
            active      false;
        }
    }

    residualAlpha   1e-16;
}

water
{
    type            pureIsothermalPhaseModel;

    diameterModel   constant;

    constantCoeffs
    {
        d           1e-4;
    }

    residualAlpha   1e-16;
}

blending
{
    default
    {
        type            linear;
        minFullyContinuousAlpha.air 0.7;
        minPartlyContinuousAlpha.air 0.3;
        minFullyContinuousAlpha.water 0.7;
        minPartlyContinuousAlpha.water 0.3;
    }

    drag
    {
        type            linear;
        minFullyContinuousAlpha.air 0.7;
        minPartlyContinuousAlpha.air 0.5;
        minFullyContinuousAlpha.water 0.7;
        minPartlyContinuousAlpha.water 0.5;
    }

    polyDrag
    {
        $drag;
    }
}

surfaceTension
{
    air_water
    {
        type    constant;
        sigma   0.07197;
    }
}

interfaceCompression
{}

drag
{
    air_dispersedIn_water
    {
        type            LogMoM;
        mode            quadrature;
        GaussHermite    10;

        drag
        {
            type        TomiyamaCorrelated;
            A           24.0;
            residualRe  1e-3;
        }
    }

    water_dispersedIn_air
    {
        type        SchillerNaumann;
        residualRe  1e-3;
    }

    air_segregatedWith_water
    {
        type    segregated;
        m       0.5;
        n       8;
    }
}

virtualMass
{
    air_dispersedIn_water
    {
        type    constantCoefficient;
        Cvm     0.5;
    }

    water_dispersedIn_air
    {
        type    constantCoefficient;
        Cvm     0.5;
    }
}

heatTransfer
{}

phaseTransfer
{}

lift
{
    air_dispersedIn_water
    {
        type            LogMoM;
        mode            quadrature;
        GaussHermite    10;

        lift
        {
            type    hyperbolic;
            C1      0.288;
            C2     -0.288;
            d       0.005;
            C       2.0;
        }
    }

    water_dispersedIn_air
    {
        type    none;
    }
}

wallLubrication
{}

turbulentDispersion
{
    air_dispersedIn_water
    {
        type            LogMoM;
        mode            quadrature;
        GaussHermite    10;

        turbulentDispersion
        {
            type            LopezDeBertodano;
            Ctd             1;
            residualAlpha   1e-16;
        }
    }

    water_dispersedIn_air
    {
        type    none;
    }
}

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
        dMax            1e-1;
        dMin            1e-4;

        p               3;
        q               2;

        continuousPhase water;

        coalescence
        {
            active                  true;
            efficiencyType          PrinceBlanch;
            frequencyType           PrinceBlanch;

            PrinceBlanchCoaEffCoeffs
            {
                sigma               0.07197;
            }

            PrinceBlanchCoaFreqCoeffs
            {
                sigma                   0.07197;
                turbulentCoalescence    true;
                buoyantCoalescence      true;
                laminarCoalescence      false;
            }

            GaussHermite            5;
        }

        breakup
        {
            active          true;
            type            LuoSvendsen;

            LuoSvendsenBreakCoeffs
            {
                sigma       0.07197;
            }

            GaussHermite    5;
            GaussLegendre   10;
        }
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

aspectRatio
{
    air_dispersedIn_water
    {
        type Wellek;
    }
}

drag
{
    air_dispersedIn_water
    {
        type            LogMoM;
        mode            quadrature;
        GaussHermite    5;

        drag
        {
            type        IshiiZuber;
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
        GaussHermite    5;

        lift
        {
            type    wallDamped;

            lift
            {
                type    Tomiyama;

                aspectRatio
                {
                    type    Wellek;
                }
            }

            wallDamping
            {
                type                linear;
                Cd                  1;
                zeroInNearWallCells true;
            }
        }
    }

    water_dispersedIn_air
    {
        type    none;
    }
}

wallLubrication
{}

wallLubrication
{
    air_dispersedIn_water
    {
        type            LogMoM;
        mode            quadrature;
        GaussHermite    5;

        wallLubrication
        {
            type    Frank;
            Cwc	    10;
            Cwd     6.8;
            p       1.7;
        }
    }

    water_dispersedIn_air
    {
        type    none;
    }
}

turbulentDispersion
{
    air_dispersedIn_water
    {
        type            LogMoM;
        mode            quadrature;
        GaussHermite    5;

        turbulentDispersion
        {
            type    Burns;
            sigma   1;
        }
    }

    water_dispersedIn_air
    {
        type    none;
    }
}

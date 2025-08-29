FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvConstraints;
}

limitp
{
    type	limitPressure;
    min		1e4;
}

limitUWater
{
    type        limitMag;
    max         10.0;
    cellZone    all;
    field       U.water;
}

limitUAir
{
    type        limitMag;
    max         10.0;
    cellZone    all;
    field       U.steam;
}

limitTSteam
{
    type        limitTemperature;
    min         VARTMIN;
    max         VARTMAX;
    cellZone    all;
    phase       steam;
}

limitTWater
{
    type        limitTemperature;
    min         VARTMIN;
    max         VARTMAX;
    cellZone    all;
    phase       water;
}

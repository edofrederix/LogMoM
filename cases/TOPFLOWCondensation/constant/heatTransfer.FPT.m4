FoamFile
{
    format      ascii;
    class       dictionary;
    object      heatTransfer;
}

steam_dispersedIn_water_inThe_water
{
    type            VARHEATMODEL;

    #include        "saturationTemperature"
    mDot            phaseChange:mDot;
    firstPhase      steam;
}

steam_dispersedIn_water_inThe_steam
{
    type            VARHEATMODEL;

    #include        "saturationTemperature"
    mDot            phaseChange:mDot;
    firstPhase      steam;
}

water_dispersedIn_steam_inThe_water
{
    type            spherical;
}

water_dispersedIn_steam_inThe_steam
{
    type            spherical;
}

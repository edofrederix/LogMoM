FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties.steam;
}

simulationType  RAS;

RAS
{
    model       VARTURBMODEL_G;

    turbulence  on;
    printCoeffs on;
}

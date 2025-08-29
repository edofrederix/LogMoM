FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties.water;
}

simulationType  RAS;

RAS
{
    model       VARTURBMODEL_L;

    turbulence  on;
    printCoeffs on;
}

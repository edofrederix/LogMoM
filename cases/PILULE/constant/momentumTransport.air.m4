FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties.air;
}

simulationType  RAS;

RAS
{
    model       VARTURBMODELAIR;

    turbulence  on;
    printCoeffs on;
}

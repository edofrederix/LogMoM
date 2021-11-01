FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application     multiphaseEulerFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         VART;

deltaT          1e-5;

writeControl    adjustableRunTime;

writeInterval   1;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   10;

runTimeModifiable yes;

adjustTimeStep  yes;

maxCo           0.5;

maxDeltaT       1;

libs
(
    "libLogMoM.so"
);

functions
{
    radialValuesOutlet
    {
        type                sets;
        libs                ("libsampling.so");
        writeControl        writeTime;
        interpolationScheme cellPoint;
        setFormat           csv;
        name                outlet;
        sets
        (
            outlet
            {
                type    lineUniform;
                axis    x;
                start   (0 VARL 0);
                end     (0.0975 VARL 0);
                nPoints 200;
            }
        );
        fields
        (
            #include "sampleFields"
        );
    }
}

FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application     foamRun;

solver          multiphaseEuler;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         20;

deltaT          1e-5;

writeControl    adjustableRunTime;

writeInterval   1;

purgeWrite      0;

writeFormat     binary;

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
    "libfieldFunctionObjects.so"
    "libmultiphaseSixDoFRigidBodyMotion.so"
);

functions
{
    #includeFunc writeObjects(d.air)
    #includeFunc writeObjects(sigma.air)

    averages
    {
        type            fieldAverage;

        writeControl    writeTime;

        log             on;

        periodicRestart on;
        restartPeriod   10.0;

        fields
        (
            #include "averagingFields"
        );
    }

    forceAir
    {
        type            movingForces;
        libs            ("libforces.so");
        log             on;
        writeControl    timeStep;
        writeInterval   1;
        patches         (wall_cylinder);
        CofR            (0 0 0);
        p               p;
        phase           air;
    }

    forceWater
    {
        type            movingForces;
        libs            ("libforces.so");
        log             on;
        writeControl    timeStep;
        writeInterval   1;
        patches         (wall_cylinder);
        CofR            (0 0 0);
        p               p;
        phase           water;
    }

    CFL
    {
        type            CourantNo;
        writeControl    writeTime;
        log             on;
    }

    outflow
    {
        type            surfaceFieldValue;
        writeFields     false;
        select          patch;
        patch           outlet;
        operation       sum;
        fields
        (
            alphaPhi.air
            alphaRhoPhi.air
            alphaPhi.water
            alphaRhoPhi.water
        );
    }

    inflow
    {
        type            surfaceFieldValue;
        writeFields     false;
        select          patch;
        patch           inlet;
        operation       sum;
        fields
        (
            alphaPhi.air
            alphaRhoPhi.air
            alphaPhi.water
            alphaRhoPhi.water
        );
    }

    #include "sample"
    #include "surface"
    #include "sixDoFRigidBodyState"
}

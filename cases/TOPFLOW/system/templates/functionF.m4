VARFUNCTIONNAME
{
    type            surfaceFieldValue;
    regionType      patch;
    name            VARPATCHNAME;
    fields          (VARFFIELD);
    operation       sum;
    weightField     VARALPHAFLUX;
    writeControl    writeTime;
    writeFields     false;
}

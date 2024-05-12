VARFUNCTIONNAME
{
    type            surfaceFieldValue;
    regionType      patch;
    name            VARPATCHNAME;
    fields          (VARFLUX);
    operation       sum;
    writeControl    writeTime;
    writeFields     false;
}

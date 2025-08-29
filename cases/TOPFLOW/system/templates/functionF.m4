VARFUNCTIONNAME
{
    type            surfaceFieldValue;
    patch           VARPATCHNAME;
    fields          (VARFFIELD);
    operation       sum;
    weightField     VARALPHAFLUX;
    writeControl    writeTime;
    writeFields     false;
}

VARFUNCTIONNAME
{
    type            surfaceFieldValue;
    patch           VARPATCHNAME;
    fields          (VARFLUX);
    operation       sum;
    writeControl    writeTime;
    writeFields     false;
}

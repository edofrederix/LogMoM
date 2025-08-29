FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      heatTransfer;
}

steam_dispersedIn_water_inThe_water
{
    type 		    LogMoM;
    mode		    quadrature;
    GaussHermite	3;

    heatTransfer
    {
        type    VARNUCORR;
        Nu      VARNU;
        Nu0     VARNU;
        d0      VARD0;
    }
}

steam_dispersedIn_water_inThe_steam
{
    type    spherical;
}

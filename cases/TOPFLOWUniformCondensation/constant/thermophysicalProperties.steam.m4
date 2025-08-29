FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      thermophysicalProperties.steam;
}

thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        nMoles      1;
        molWeight   18.0153;
    }
    thermodynamics
    {
        Cp          VARCPVAP;
        Hf          0;
        Tref        VARTR;
        Hsref       VARHVAP;
    }
    equationOfState
    {
        rho         VARRHOVAP;
    }
    transport
    {
        mu          VARMUVAP;
        Pr          VARPRVAP;
    }
}

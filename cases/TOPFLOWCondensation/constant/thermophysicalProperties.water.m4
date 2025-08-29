FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      thermophysicalProperties.water;
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
        molWeight   18.0153;
    }
    equationOfState
    {
        rho         VARRHOLIQ;
    }
    thermodynamics
    {
        Cp          VARCPLIQ;
        Hf          0;
        Tref        VARTR;
        Hsref       VARHLIQ;
    }
    transport
    {
        mu          VARMULIQ;
        Pr          VARPRLIQ;
    }
}

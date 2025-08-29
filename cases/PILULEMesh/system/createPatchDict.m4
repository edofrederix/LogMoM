FoamFile
{
    format      ascii;
    class       dictionary;
    object      createPatchDict;
}

patches
(
    {
        name        VARPATCHNAME;

        patchInfo
        {
            type    patch;
        }

        constructFrom set;
        set         symm;
    }
);

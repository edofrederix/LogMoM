VARPHASEPAIR
{
    type    wallDamped;

    lift
    {
        type    Tomiyama;
        Cl      0.288;

        aspectRatio
        {
            type    constant;
            E0      1;
        }
    }

    wallDamping
    {
        type                linear;
        Cd                  1;
        zeroInNearWallCells true;
    }
}

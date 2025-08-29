FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      mirrorMeshDict;
}

planeType       pointAndNormal;

point           (0 0 0);
normal          (-1 0 0);

planeTolerance  1e-06;

VARPHASENAME
{
    type            pureIsothermalPhaseModel;

    diameterModel   velocityGroup;

    velocityGroupCoeffs
    {
        populationBalance    air;

        shapeModel           spherical;

        sizeGroups
        (
            @include "sizeGroups.VARPHASENAME"
        );
    }

    residualAlpha   1e-16;
}

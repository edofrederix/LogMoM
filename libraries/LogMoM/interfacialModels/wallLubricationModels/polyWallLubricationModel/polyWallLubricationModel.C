#include "polyWallLubricationModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyWallLubricationModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(polyWallLubricationModel, 0);
    defineRunTimeSelectionTable(polyWallLubricationModel, dictionary);
}

const Foam::dimensionSet Foam::polyWallLubricationModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyWallLubricationModel::polyWallLubricationModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    wallDependentModel(interface.mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyWallLubricationModel::~polyWallLubricationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::blendedPolyWallLubricationModel::F
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyWallLubricationModel::F,
            "F",
            polyWallLubricationModel::dimF,
            true,
            gamma
        );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::blendedPolyWallLubricationModel::Ff
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyWallLubricationModel::Ff,
            "Ff",
            polyWallLubricationModel::dimF*dimArea,
            true,
            gamma
        );
}


// ************************************************************************* //

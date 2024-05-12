#include "defaultPolyWallLubricationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(defaultPolyWallLubricationModel, 0);
    addToRunTimeSelectionTable
    (
        polyWallLubricationModel,
        defaultPolyWallLubricationModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::defaultPolyWallLubricationModel::
defaultPolyWallLubricationModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    polyWallLubricationModel(dict, interface),
    wallLubricationModelPtr_
    (
        wallLubricationModel::New(dict, interface, false).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::defaultPolyWallLubricationModel::
~defaultPolyWallLubricationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::defaultPolyWallLubricationModel::F
(
    const label gamma
) const
{
    return wallLubricationModelPtr_->F();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::wallLubricationModels::defaultPolyWallLubricationModel::Ff
(
    const label gamma
) const
{
    return wallLubricationModelPtr_->Ff();
}

// ************************************************************************* //

#include "LogMoMWallLubricationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(LogMoMWallLubricationModel, 0);
    addToRunTimeSelectionTable(wallLubricationModel, LogMoMWallLubricationModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::LogMoMWallLubricationModel::LogMoMWallLubricationModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    wallLubricationModel(dict, interface),
    wallLubricationModelPtr_
    (
        wallLubricationModel::New
        (
            dict.subDict("wallLubrication"),
            interface,
            false
        ).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::LogMoMWallLubricationModel::
~LogMoMWallLubricationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::LogMoMWallLubricationModel::F() const
{
    return wallLubricationModelPtr_->F();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::wallLubricationModels::LogMoMWallLubricationModel::Ff()
const
{
    return wallLubricationModelPtr_->Ff();
}

// ************************************************************************* //

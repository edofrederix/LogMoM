#include "LogMoMWallLubricationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedWallLubricationModel.H"
#include "fvcFlux.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(LogMoMWallLubricationModel, 0);
    addToRunTimeSelectionTable
    (
        wallLubricationModel,
        LogMoMWallLubricationModel,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::LogMoMWallLubricationModel::FiModel() const
{
    return
        refCast<const dispersedWallLubricationModel>
        (
            wallLubricationModel_()
        ).Fi();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::LogMoMWallLubricationModel::
LogMoMWallLubricationModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedWallLubricationModel(dict, interface),
    LogMoMInterfacialModel(dict, interface),
    wallLubricationModel_
    (
        wallLubricationModel::New
        (
            dict.subDict("wallLubrication"),
            interface,
            false
        )
    )
{
    if (!isA<dispersedWallLubricationModel>(wallLubricationModel_()))
    {
        FatalErrorInFunction
            << "The sub-wall-lubrication-model of a " << type()
            << " wall lubrication model must be for a dispersed configuration"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::LogMoMWallLubricationModel::
~LogMoMWallLubricationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::LogMoMWallLubricationModel::Fi() const
{
    return this->evaluate
    (
        gamma(),
        &LogMoMWallLubricationModel::FiModel,
        *this
    );
}

// ************************************************************************* //

#include "LogMoMPolyWallLubricationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedWallLubricationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(LogMoMPolyWallLubricationModel, 0);
    addToRunTimeSelectionTable
    (
        polyWallLubricationModel,
        LogMoMPolyWallLubricationModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::LogMoMPolyWallLubricationModel::LogMoMPolyWallLubricationModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    polyWallLubricationModel(dict, interface),
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

Foam::wallLubricationModels::LogMoMPolyWallLubricationModel::~LogMoMPolyWallLubricationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModels::LogMoMPolyWallLubricationModel::Fim() const
{
    return
        refCast<const dispersedWallLubricationModel>
        (
            wallLubricationModel_()
        ).Fi();
}

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModels::LogMoMPolyWallLubricationModel::Fi() const
{
    return
        refCast<const dispersedWallLubricationModel>
        (
            wallLubricationModel_()
        ).Fi();
}

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::LogMoMPolyWallLubricationModel::Fi
(
    const label gamma
) const
{
    return this->evaluate(gamma, &LogMoMPolyWallLubricationModel::Fim, *this);
}

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModels::LogMoMPolyWallLubricationModel::F
(
    const label gamma
) const
{
    return interface_.dispersed()*Fi(gamma);
}

Foam::tmp<Foam::surfaceScalarField> Foam::wallLubricationModels::LogMoMPolyWallLubricationModel::Ff
(
    const label gamma
) const
{
    return wallLubricationModel_->Ff();
}

// ************************************************************************* //

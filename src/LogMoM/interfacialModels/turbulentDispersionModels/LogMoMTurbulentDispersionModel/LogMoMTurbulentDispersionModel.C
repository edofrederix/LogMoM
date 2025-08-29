#include "LogMoMTurbulentDispersionModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedTurbulentDispersionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(LogMoMTurbulentDispersionModel, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        LogMoMTurbulentDispersionModel,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::LogMoMTurbulentDispersionModel::DModel() const
{
    return
        refCast<const dispersedTurbulentDispersionModel>
        (
            turbulentDispersionModel_()
        ).D();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::LogMoMTurbulentDispersionModel::
LogMoMTurbulentDispersionModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedTurbulentDispersionModel(dict, interface),
    LogMoMInterfacialModel(dict, interface),
    turbulentDispersionModel_
    (
        turbulentDispersionModel::New
        (
            dict.subDict("turbulentDispersion"),
            interface,
            false
        )
    )
{
    if (!isA<dispersedTurbulentDispersionModel>(turbulentDispersionModel_()))
    {
        FatalErrorInFunction
            << "The sub-turbulent-dispersion-model of a " << type()
            << " turbulent dispersion model must be for a dispersed "
            << " configuration"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::LogMoMTurbulentDispersionModel
::~LogMoMTurbulentDispersionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::LogMoMTurbulentDispersionModel::D() const
{
    return this->evaluate
    (
        gamma(),
        &LogMoMTurbulentDispersionModel::DModel,
        *this
    );
}

// ************************************************************************* //

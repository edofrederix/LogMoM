#include "LogMoMPolyTurbulentDispersionModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedTurbulentDispersionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(LogMoMPolyTurbulentDispersionModel, 0);
    addToRunTimeSelectionTable
    (
        polyTurbulentDispersionModel,
        LogMoMPolyTurbulentDispersionModel,
        dictionary
    );
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::LogMoMPolyTurbulentDispersionModel::
LogMoMPolyTurbulentDispersionModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    polyTurbulentDispersionModel(dict, interface),
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

Foam::turbulentDispersionModels::LogMoMPolyTurbulentDispersionModel
::~LogMoMPolyTurbulentDispersionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::LogMoMPolyTurbulentDispersionModel::D() const
{
    return
        refCast<const dispersedTurbulentDispersionModel>
        (
            turbulentDispersionModel_()
        ).D();
}

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::LogMoMPolyTurbulentDispersionModel::D
(
    const label gamma
) const
{
    return this->evaluate(gamma, &LogMoMPolyTurbulentDispersionModel::D, *this);
}

// ************************************************************************* //

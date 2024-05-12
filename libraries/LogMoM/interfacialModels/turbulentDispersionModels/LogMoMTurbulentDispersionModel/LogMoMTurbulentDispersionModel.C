#include "LogMoMTurbulentDispersionModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedDragModel.H"

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::LogMoMTurbulentDispersionModel::LogMoMTurbulentDispersionModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    turbulentDispersionModel(dict, interface),
    turbulentDispersionModelPtr_
    (
        turbulentDispersionModel::New
        (
            dict.subDict("turbulentDispersion"),
            interface,
            false
        ).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::LogMoMTurbulentDispersionModel::~LogMoMTurbulentDispersionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::LogMoMTurbulentDispersionModel::D() const
{
    return turbulentDispersionModelPtr_->D();
}


// ************************************************************************* //

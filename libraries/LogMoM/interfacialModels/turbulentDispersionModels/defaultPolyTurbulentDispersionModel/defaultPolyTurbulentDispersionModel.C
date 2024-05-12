#include "defaultPolyTurbulentDispersionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(defaultPolyTurbulentDispersionModel, 0);
    addToRunTimeSelectionTable
    (
        polyTurbulentDispersionModel,
        defaultPolyTurbulentDispersionModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::defaultPolyTurbulentDispersionModel::
defaultPolyTurbulentDispersionModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    polyTurbulentDispersionModel(dict, interface),
    turbulentDispersionModelPtr_
    (
        turbulentDispersionModel::New(dict, interface, false).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::defaultPolyTurbulentDispersionModel::
~defaultPolyTurbulentDispersionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::defaultPolyTurbulentDispersionModel::D
(
    const label gamma
) const
{
    return turbulentDispersionModelPtr_->D();
}


// ************************************************************************* //

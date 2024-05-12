#include "defaultPolyVirtualMassModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(defaultPolyVirtualMassModel, 0);
    addToRunTimeSelectionTable
    (
        polyVirtualMassModel,
        defaultPolyVirtualMassModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::defaultPolyVirtualMassModel::
defaultPolyVirtualMassModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    polyVirtualMassModel(dict, interface, registerObject),
    virtualMassModelPtr_
    (
        virtualMassModel::New(dict, interface, false, registerObject).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModels::defaultPolyVirtualMassModel::
~defaultPolyVirtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::defaultPolyVirtualMassModel::K
(
    const label gamma
) const
{
    return virtualMassModelPtr_->K();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::virtualMassModels::defaultPolyVirtualMassModel::Kf
(
    const label gamma
) const
{
    return virtualMassModelPtr_->Kf();
}

// ************************************************************************* //

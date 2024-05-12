#include "LogMoMVirtualMassModel.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(LogMoMVirtualMassModel, 0);
    addToRunTimeSelectionTable
    (
        virtualMassModel,
        LogMoMVirtualMassModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::LogMoMVirtualMassModel::LogMoMVirtualMassModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    virtualMassModel(dict, interface, false),
    virtualMassModelPtr_
    (
        virtualMassModel::New
        (
            dict.subDict("virtualMass"),
            interface,
            false,
            registerObject
        ).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModels::LogMoMVirtualMassModel::~LogMoMVirtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::LogMoMVirtualMassModel::K() const
{
    return virtualMassModelPtr_->K();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::virtualMassModels::LogMoMVirtualMassModel::Kf() const
{
    return virtualMassModelPtr_->Kf();
}

// ************************************************************************* //

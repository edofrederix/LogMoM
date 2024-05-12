#include "LogMoMPolyVirtualMassModel.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedVirtualMassModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(LogMoMPolyVirtualMassModel, 0);
    addToRunTimeSelectionTable
    (
        polyVirtualMassModel,
        LogMoMPolyVirtualMassModel,
        dictionary
    );
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::LogMoMPolyVirtualMassModel::LogMoMPolyVirtualMassModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    polyVirtualMassModel(dict, interface, registerObject),
    LogMoMInterfacialModel(dict, interface),
    virtualMassModelPtr_
    (
        virtualMassModel::New
        (
            dict.subDict("virtualMass"),
            interface,
            false,
            registerObject
        )
    )
{
    if (!isA<dispersedVirtualMassModel>(virtualMassModelPtr_()))
    {
        FatalErrorInFunction
            << "The sub-turbulent-dispersion-model of a " << type()
            << " turbulent dispersion model must be for a dispersed "
            << " configuration"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModels::LogMoMPolyVirtualMassModel::
~LogMoMPolyVirtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::LogMoMPolyVirtualMassModel::Cvm() const
{
    return refCast<const dispersedVirtualMassModel>
    (
        virtualMassModelPtr_()
    ).Cvm();
}

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::LogMoMPolyVirtualMassModel::Ki() const
{
    return refCast<const dispersedVirtualMassModel>
    (
        virtualMassModelPtr_()
    ).Ki();
}

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::LogMoMPolyVirtualMassModel::Ki
(
    const label gamma
) const
{
    return
        this->evaluate(gamma, &LogMoMPolyVirtualMassModel::Cvm, *this)
      * interface_.continuous().rho();
}

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::LogMoMPolyVirtualMassModel::K
(
    const label gamma
) const
{
    return interface_.dispersed()*Ki(gamma);
}

Foam::tmp<Foam::surfaceScalarField>
Foam::virtualMassModels::LogMoMPolyVirtualMassModel::Kf
(
    const label gamma
) const
{
    return virtualMassModelPtr_->Kf();
}

// ************************************************************************* //

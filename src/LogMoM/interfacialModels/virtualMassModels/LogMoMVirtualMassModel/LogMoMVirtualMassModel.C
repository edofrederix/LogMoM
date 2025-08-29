#include "LogMoMVirtualMassModel.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"
#include "dispersedVirtualMassModel.H"

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
    dispersedVirtualMassModel(dict, interface, registerObject),
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

Foam::virtualMassModels::LogMoMVirtualMassModel::~LogMoMVirtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::LogMoMVirtualMassModel::Cvm() const
{
    return refCast<const dispersedVirtualMassModel>
    (
        virtualMassModelPtr_()
    ).Cvm();
}

Foam::tmp<Foam::volScalarField>
Foam::virtualMassModels::LogMoMVirtualMassModel::Ki() const
{
    return
        this->evaluate(gamma(), &LogMoMVirtualMassModel::Cvm, *this)
      * interface_.continuous().rho();
}

// ************************************************************************* //

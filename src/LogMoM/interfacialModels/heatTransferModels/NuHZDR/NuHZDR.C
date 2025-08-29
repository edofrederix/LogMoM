#include "NuHZDR.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
// #include "interfaceSaturationTemperatureModel.H"
// #include "phaseInterfaceKey.H"
// #include "HeatTransferPhaseSystem.H"
// #include "MomentumTransferPhaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(NuHZDR, 0);
    addToRunTimeSelectionTable(heatTransferModel, NuHZDR, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::NuHZDR::NuHZDR
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    heatTransferModel(dict,interface,registerObject),
    interface_
    (
        interface.modelCast<heatTransferModel, dispersedPhaseInterface>()
    ),
    saturationModelPtr_
    (
        saturationTemperatureModel::New
        (
            "saturationTemperature",
            dict
        ).ptr()
    ),
    mDotName_(dict.lookup("mDot")),
    firstPhaseName_(dict.lookup("firstPhase"))
{
    if
    (
        firstPhaseName_ != interface_.dispersed().name()
     && firstPhaseName_ != interface_.continuous().name()
    )
    {
        FatalErrorInFunction
            << "Invalid first phase. Available phases are "
            << interface_.dispersed().name() << " and "
            << interface_.continuous().name() << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::NuHZDR::~NuHZDR()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::phaseCompressible::momentumTransportModel&
Foam::heatTransferModels::NuHZDR::continuousTurbulence() const
{
    return
        interface_.continuous().mesh()
       .lookupObject<phaseCompressible::momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                interface_.continuous().name()
            )
        );
}

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::NuHZDR::K(const scalar residualAlpha) const
{
    const scalar pi(constant::mathematical::pi);

    const phaseModel& dispPhase = interface_.dispersed();
    const phaseModel& contPhase = interface_.continuous();

    const fvMesh& mesh = dispPhase.mesh();

    const volScalarField& T = contPhase.thermo().T();

    const volScalarField& rhod = dispPhase.rho();
    const volScalarField& rhoc = contPhase.rho();

    const volScalarField& p = dispPhase.fluidThermo().p();

    const volScalarField Tsat(saturationModelPtr_->Tsat(p));

    const volScalarField hafc(contPhase.fluidThermo().ha(p, Tsat));
    const volScalarField hafd(dispPhase.fluidThermo().ha(p, Tsat));
    const volScalarField hac(contPhase.fluidThermo().ha());
    const volScalarField had(dispPhase.fluidThermo().ha());

    // The mDot field is redefined here to be positive for net transfer from
    // dispersed to continuous. This involves a sign change if the first phase
    // of the phase change model is not the dispersed phase.

    const volScalarField mDot
    (
        IOobject("mDot", mesh.time().name(), mesh),
        (dispPhase.name() == firstPhaseName_ ? 1.0 : -1.0)
      * mesh.lookupObject<volScalarField::Internal>(mDotName_),
        calculatedFvPatchScalarField::typeName
    );

    const volScalarField L
    (
        neg0(mDot)*hafc + pos(mDot)*hac
      - pos0(mDot)*hafd - neg(mDot)*had
    );

    const volScalarField Re(interface_.Re());
    const volScalarField Pr(interface_.Pr());
    const volScalarField k(continuousTurbulence().k());
    const volScalarField epsilon(continuousTurbulence().epsilon());

    const volScalarField lt(pow(0.09, 3.0/4.0)*pow(k, 3.0/2.0)/epsilon);
    const volScalarField ut(pow(0.09, 1.0/4.0)*sqrt(k));

    const volScalarField muc(contPhase.fluidThermo().mu());
    const volScalarField Cpc(contPhase.thermo().Cp());

    const volScalarField Ret(rhoc*lt*ut/muc);
    const volScalarField Pe(Re*Pr);
    const volScalarField Pet(Ret*Pr);
    const volScalarField Ja(rhoc*Cpc*mag(Tsat - T)/(rhod*L));

    const volScalarField d(interface_.dispersed().d());

    const volScalarField Nu
    (
        12.0/pi*Ja
      + 2.0*sqrt(Pe/pi)
      + 2.0*sqrt(Pet/pi)*d/lt
    );

    return
        6.0*max(interface_.dispersed(), residualAlpha)
      * interface_.continuous().thermo().kappa()
      * Nu
      / sqr(d);
}


// ************************************************************************* //
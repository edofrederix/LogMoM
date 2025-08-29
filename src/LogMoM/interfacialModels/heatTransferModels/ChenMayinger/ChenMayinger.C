#include "ChenMayinger.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(ChenMayinger, 0);
    addToRunTimeSelectionTable(heatTransferModel, ChenMayinger, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::ChenMayinger::ChenMayinger
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::ChenMayinger::~ChenMayinger()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::ChenMayinger::K(const scalar residualAlpha) const
{
    volScalarField Nu(0.185*pow(interface_.Re(), 0.7)*sqrt(interface_.Pr()));

    return
        6
        *max(interface_.dispersed(), residualAlpha)
        *interface_.continuous().thermo().kappa()
        *Nu
        /sqr(interface_.dispersed().d());
}


// ************************************************************************* //
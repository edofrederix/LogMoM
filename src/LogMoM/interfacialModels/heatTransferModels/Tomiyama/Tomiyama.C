#include "Tomiyama.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(Tomiyama, 0);
    addToRunTimeSelectionTable(heatTransferModel, Tomiyama, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::Tomiyama::Tomiyama
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

Foam::heatTransferModels::Tomiyama::~Tomiyama()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::Tomiyama::K(const scalar residualAlpha) const
{
    volScalarField Nu(2 + 0.15*pow(interface_.Re(), 0.8)*sqrt(interface_.Pr()));

    return
        6
        *max(interface_.dispersed(), residualAlpha)
        *interface_.continuous().thermo().kappa()
        *Nu
        /sqr(interface_.dispersed().d());
}


// ************************************************************************* //
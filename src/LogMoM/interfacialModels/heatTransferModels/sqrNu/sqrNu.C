#include "sqrNu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(sqrNu, 0);
    addToRunTimeSelectionTable(heatTransferModel, sqrNu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::sqrNu::sqrNu
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
    Nu0_("Nu0", dimless, dict),
    d0_("d0", dimLength, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::sqrNu::~sqrNu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::sqrNu::K(const scalar residualAlpha) const
{
    return
        6.0
      * max(interface_.dispersed(), residualAlpha)
      * interface_.continuous().thermo().kappa()
      * Nu0_
      / sqr(d0_);
}


// ************************************************************************* //
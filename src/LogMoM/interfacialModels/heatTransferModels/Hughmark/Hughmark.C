#include "Hughmark.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(Hughmark, 0);
    addToRunTimeSelectionTable(heatTransferModel, Hughmark, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::Hughmark::Hughmark
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

Foam::heatTransferModels::Hughmark::~Hughmark()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::Hughmark::K(const scalar residualAlpha) const
{
    volScalarField Rep(interface_.Re());
    volScalarField Nu(2 + 0.6*sqrt(Rep)*cbrt(interface_.Pr()));

    const fvMesh& mesh_(interface_.dispersed().mesh());
    volScalarField Recrit
    (
        IOobject
        (
             "Recrit",
             mesh_.time().name(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 776.06)
    );

    if (Rep > Recrit)
    {
        Nu = 2 + 0.27*pow(Rep, 0.62)*cbrt(interface_.Pr());
    }

    return
        6
        *max(interface_.dispersed(), residualAlpha)
        *interface_.continuous().thermo().kappa()
        *Nu
        /sqr(interface_.dispersed().d());
}


// ************************************************************************* //
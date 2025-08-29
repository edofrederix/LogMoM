#include "LogMoMHeatTransferModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(LogMoMHeatTransferModel, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel,
        LogMoMHeatTransferModel,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::LogMoMHeatTransferModel::KdSqrPow() const
{
    return
        refCast<const heatTransferModel>
        (
            heatTransferModelPtr_()
        ).K(residualAlpha_)
      * pow(interface_.dispersed().d(), 2.0 - pow_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::LogMoMHeatTransferModel::
LogMoMHeatTransferModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    heatTransferModel(dict, interface, registerObject),
    LogMoMInterfacialModel(dict, interface),
    interface_
    (
        interface.modelCast<heatTransferModel, dispersedPhaseInterface>()
    ),
    heatTransferModelPtr_
    (
        heatTransferModel::New
        (
            dict.subDict("heatTransfer"),
            interface,
            false
        ).ptr()
    ),
    pow_(dict.lookupOrDefault<scalar>("power", 0.5))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::LogMoMHeatTransferModel::
~LogMoMHeatTransferModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::LogMoMHeatTransferModel::K
(
    const scalar residualAlpha
) const
{
    const_cast<scalar&>(residualAlpha_) = residualAlpha;

    // We need to take the gamma'th moment of K which is proportional to Nu/d^2
    // and Nu is proportional to d^pow, so K ~ d^(pow-2). Instead, we take the
    // (gamma+2-pow)'th moment of the scaled function K*d^(2-pow), so that the
    // evaluation  is more accurate.

    tmp<volScalarField> tK =
        this->evaluate
        (
            gamma() + pow_ - 2.0,
            gamma(),
            &LogMoMHeatTransferModel::KdSqrPow,
            *this
        );

    return tK;
}

// ************************************************************************* //

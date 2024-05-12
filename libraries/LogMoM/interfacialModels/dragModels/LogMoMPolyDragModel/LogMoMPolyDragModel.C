#include "LogMoMPolyDragModel.H"
#include "addToRunTimeSelectionTable.H"
#include "swarmCorrection.H"
#include "threeMomentLogNormal.H"
#include "noSwarm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(LogMoMPolyDragModel, 0);
    addToRunTimeSelectionTable(polyDragModel, LogMoMPolyDragModel, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::LogMoMPolyDragModel::LogMoMPolyDragModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    polyDragModel(dict, interface, registerObject),
    LogMoMInterfacialModel(dict, interface),
    dragModelPtr_
    (
        dragModel::New
        (
            dict.subDict("drag"),
            interface,
            false,
            registerObject
        ).ptr()
    ),
    swarmCorrection_
    (
        dict.found("swarmCorrection")
      ? swarmCorrection::New(dict.subDict("swarmCorrection"), interface).ptr()
      : new swarmCorrections::noSwarm(dict, interface)
    )
{
    if (!isA<dispersedDragModel>(dragModelPtr_()))
    {
        FatalErrorInFunction
            << "The sub-drag-model of a " << type()
            << " drag model must be for a dispersed configuration"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::LogMoMPolyDragModel::~LogMoMPolyDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMPolyDragModel::CdRe() const
{
    return refCast<const dispersedDragModel>(dragModelPtr_()).CdRe();
}

Foam::tmp<Foam::volScalarField>
Foam::dragModels::LogMoMPolyDragModel::CdReBydSqr() const
{
    return this->CdRe()/sqr(interface_.dispersed().d());
}

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMPolyDragModel::Ki() const
{
    return refCast<const dispersedDragModel>(dragModelPtr_()).Ki();
}

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMPolyDragModel::Ki
(
    const label gamma
) const
{
    return
        0.75
      * this->evaluate(gamma, &LogMoMPolyDragModel::CdReBydSqr, *this)
      * swarmCorrection_->Cs()
      * interface_.continuous().rho()
      * interface_.continuous().thermo().nu();
}

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMPolyDragModel::K
(
    const label gamma
) const
{
    return
        max
        (
            interface_.dispersed(),
            interface_.dispersed().residualAlpha()
        )*Ki(gamma);
}

Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::LogMoMPolyDragModel::Kf
(
    const label gamma
) const
{
    return dragModelPtr_->Kf();
}

// ************************************************************************* //

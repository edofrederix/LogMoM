#include "LogMoMDragModel.H"
#include "addToRunTimeSelectionTable.H"
#include "swarmCorrection.H"
#include "noSwarm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(LogMoMDragModel, 0);
    addToRunTimeSelectionTable(dragModel, LogMoMDragModel, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::LogMoMDragModel::LogMoMDragModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
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

Foam::dragModels::LogMoMDragModel::~LogMoMDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMDragModel::CdRe() const
{
    return refCast<const dispersedDragModel>(dragModelPtr_()).CdRe();
}

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMDragModel::Ki() const
{
    // We need to take the gamma'th moment of Cd*Re/d^2. So, instead, we take
    // the (gamma - 2)'th moment of Cd*Re, so that the centered Gauss-Hermite
    // integration is more accurate.

    return
        0.75
      * this->evaluate(gamma() - 2, gamma(), &LogMoMDragModel::CdRe, *this)
      * swarmCorrection_->Cs()
      * interface_.continuous().rho()
      * interface_.continuous().fluidThermo().nu();
}

// ************************************************************************* //

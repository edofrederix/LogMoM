#include "LogMoMDragModel.H"
#include "addToRunTimeSelectionTable.H"

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
    dragModel(dict, interface, false),
    dragModelPtr_
    (
        dragModel::New
        (
            dict.subDict("drag"),
            interface,
            false,
            registerObject
        ).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::LogMoMDragModel::~LogMoMDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::LogMoMDragModel::K() const
{
    return dragModelPtr_->K();
}

Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::LogMoMDragModel::Kf()
const
{
    return dragModelPtr_->Kf();
}

// ************************************************************************* //

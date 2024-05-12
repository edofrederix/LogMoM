#include "defaultPolyDragModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(defaultPolyDragModel, 0);
    addToRunTimeSelectionTable(polyDragModel, defaultPolyDragModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::defaultPolyDragModel::defaultPolyDragModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    polyDragModel(dict, interface, registerObject),
    dragModelPtr_
    (
        dragModel::New(dict, interface, false, registerObject).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::defaultPolyDragModel::~defaultPolyDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::defaultPolyDragModel::K
(
    const label gamma
) const
{
    return dragModelPtr_->K();
}

Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::defaultPolyDragModel::Kf
(
    const label gamma
) const
{
    return dragModelPtr_->Kf();
}

// ************************************************************************* //

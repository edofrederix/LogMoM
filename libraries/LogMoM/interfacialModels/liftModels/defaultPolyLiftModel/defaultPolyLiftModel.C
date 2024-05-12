#include "defaultPolyLiftModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(defaultPolyLiftModel, 0);
    addToRunTimeSelectionTable(polyLiftModel, defaultPolyLiftModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::defaultPolyLiftModel::defaultPolyLiftModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    polyLiftModel(dict, interface),
    liftModelPtr_
    (
        liftModel::New(dict, interface, false).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::defaultPolyLiftModel::~defaultPolyLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::liftModels::defaultPolyLiftModel::F
(
    const label gamma
) const
{
    return liftModelPtr_->F();
}

Foam::tmp<Foam::surfaceScalarField> Foam::liftModels::defaultPolyLiftModel::Ff
(
    const label gamma
) const
{
    return liftModelPtr_->Ff();
}

// ************************************************************************* //

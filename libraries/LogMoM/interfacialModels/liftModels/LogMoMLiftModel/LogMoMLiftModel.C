#include "LogMoMLiftModel.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(LogMoMLiftModel, 0);
    addToRunTimeSelectionTable(liftModel, LogMoMLiftModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::LogMoMLiftModel::LogMoMLiftModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    liftModel(dict, interface),
    liftModelPtr_
    (
        liftModel::New(dict.subDict("lift"), interface, false).ptr()
    )
{
    if (dict.found("wallDamping"))
        wallDampingModelPtr_ =
            wallDampingModel::New(dict.subDict("wallDamping"), interface);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::LogMoMLiftModel::~LogMoMLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::liftModels::LogMoMLiftModel::F() const
{
    if (wallDampingModelPtr_.valid())
    {
        return
            liftModelPtr_->F()
          * wallDampingModelPtr_->damping();
    }
    else
    {
        return liftModelPtr_->F();
    }
}

Foam::tmp<Foam::surfaceScalarField> Foam::liftModels::LogMoMLiftModel::Ff()
const
{
    if (wallDampingModelPtr_.valid())
    {
        return
            liftModelPtr_->Ff()
          * wallDampingModelPtr_->dampingf();
    }
    else
    {
        return liftModelPtr_->Ff();
    }
}


// ************************************************************************* //

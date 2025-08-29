#include "LogMoMLiftModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcCurl.H"
#include "dispersedLiftModel.H"
#include "fvcFlux.H"

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
    dispersedLiftModel(dict, interface),
    LogMoMInterfacialModel(dict, interface),
    liftModelPtr_
    (
        liftModel::New(dict.subDict("lift"), interface, false).ptr()
    )
{
    if (!isA<dispersedLiftModel>(liftModelPtr_()))
    {
        FatalErrorInFunction
            << "The sub-lift-model of a " << type()
            << " lift model must be for a dispersed configuration"
            << exit(FatalError);
    }

    if (dict.found("wallDamping"))
        wallDampingModelPtr_ =
            wallDampingModel::New(dict.subDict("wallDamping"), interface);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::LogMoMLiftModel::~LogMoMLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::LogMoMLiftModel::Cl() const
{
    return refCast<const dispersedLiftModel>(liftModelPtr_()).Cl();
}

Foam::tmp<Foam::volVectorField> Foam::liftModels::LogMoMLiftModel::Fi() const
{
    return
        this->evaluate(gamma(), &LogMoMLiftModel::Cl, *this)
      * interface_.continuous().rho()
      * (
            interface_.Ur() ^ fvc::curl(interface_.continuous().U())
        );
}

// ************************************************************************* //

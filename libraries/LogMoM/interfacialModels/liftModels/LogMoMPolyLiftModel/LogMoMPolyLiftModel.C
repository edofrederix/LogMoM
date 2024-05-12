#include "LogMoMPolyLiftModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcCurl.H"
#include "dispersedLiftModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(LogMoMPolyLiftModel, 0);
    addToRunTimeSelectionTable(polyLiftModel, LogMoMPolyLiftModel, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::LogMoMPolyLiftModel::LogMoMPolyLiftModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    polyLiftModel(dict, interface),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::LogMoMPolyLiftModel::~LogMoMPolyLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::liftModels::LogMoMPolyLiftModel::Cl() const
{
    return refCast<const dispersedLiftModel>(liftModelPtr_()).Cl();
}

Foam::tmp<Foam::volVectorField>
Foam::liftModels::LogMoMPolyLiftModel::ClUr() const
{
    return
        refCast<const dispersedLiftModel>(liftModelPtr_()).Cl()
      * interface_.Ur();
}

Foam::tmp<Foam::volVectorField>
Foam::liftModels::LogMoMPolyLiftModel::Fi() const
{
    return refCast<const dispersedLiftModel>(liftModelPtr_()).Fi();
}

Foam::tmp<Foam::volVectorField> Foam::liftModels::LogMoMPolyLiftModel::Fi
(
    const label gamma
) const
{
    return
        this->evaluate(gamma, &LogMoMPolyLiftModel::ClUr, *this)
      * interface_.continuous().rho()
      ^ fvc::curl(interface_.continuous().U());
}

Foam::tmp<Foam::volVectorField> Foam::liftModels::LogMoMPolyLiftModel::F
(
    const label gamma
) const
{
    return interface_.dispersed()*Fi(gamma);
}

Foam::tmp<Foam::surfaceScalarField> Foam::liftModels::LogMoMPolyLiftModel::Ff
(
    const label gamma
) const
{
    return liftModelPtr_->Ff();
}

// ************************************************************************* //

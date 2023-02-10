#include "breakupModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(breakupModel, 0);
    defineRunTimeSelectionTable(breakupModel, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::breakupModel>
Foam::breakupModel::New
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
{
    word breakupModelType(dict.lookup("type"));

    Info<< "Selecting breakup model for "
        << breakupModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(breakupModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown breakup model type "
            << breakupModelType << endl << endl
            << "Valid breakup model types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(pair, dict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModel::breakupModel
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    pair_(pair.dispersed(), pair.continuous()),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::breakupModel::~breakupModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressible::momentumTransportModel&
Foam::breakupModel::continuousTurbulence() const
{
    return
        pair_.mesh().lookupObject<phaseCompressible::momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                pair_.continuous().name()
            )
        );
}

// ************************************************************************* //

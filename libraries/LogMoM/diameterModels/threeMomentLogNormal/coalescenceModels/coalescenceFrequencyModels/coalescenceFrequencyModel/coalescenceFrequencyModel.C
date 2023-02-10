#include "coalescenceFrequencyModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coalescenceFrequencyModel, 0);
    defineRunTimeSelectionTable(coalescenceFrequencyModel, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coalescenceFrequencyModel>
Foam::coalescenceFrequencyModel::New
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
{
    word coalescenceFrequencyModelType(dict.lookup("frequencyType"));

    Info<< "Selecting coalescence frequency model for "
        << coalescenceFrequencyModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(coalescenceFrequencyModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown coalescence frequency model type "
            << coalescenceFrequencyModelType << endl << endl
            << "Valid coalescence frequency model types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(pair, dict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModel::coalescenceFrequencyModel
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    pair_(pair.dispersed(), pair.continuous()),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModel::~coalescenceFrequencyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressible::momentumTransportModel&
Foam::coalescenceFrequencyModel::continuousTurbulence() const
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

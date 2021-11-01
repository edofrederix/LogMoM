#include "coalescenceEfficiencyModel.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coalescenceEfficiencyModel, 0);
    defineRunTimeSelectionTable(coalescenceEfficiencyModel, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coalescenceEfficiencyModel>
Foam::coalescenceEfficiencyModel::New
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
{
    word coalescenceEfficiencyModelType(dict.lookup("efficiencyType"));

    Info<< "Selecting coalescence efficiency model for "
        << coalescenceEfficiencyModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(coalescenceEfficiencyModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown coalescence efficiency model type "
            << coalescenceEfficiencyModelType << endl << endl
            << "Valid coalescence efficiency model types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(pair, dict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModel::coalescenceEfficiencyModel
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    pair_(pair.dispersed(), pair.continuous()),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModel::~coalescenceEfficiencyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressible::momentumTransportModel&
Foam::coalescenceEfficiencyModel::continuousTurbulence() const
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

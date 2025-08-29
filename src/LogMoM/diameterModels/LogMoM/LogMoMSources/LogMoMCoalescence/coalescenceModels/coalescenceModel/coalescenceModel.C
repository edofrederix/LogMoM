#include "coalescenceModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coalescenceModel, 0);
    defineRunTimeSelectionTable(coalescenceModel, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coalescenceModel> Foam::coalescenceModel::New
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
{
    word coalescenceModelType(dict.lookup("type"));

    Info<< "Selecting coalescence model for " << coalescenceModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(coalescenceModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown coalescence model "
            << coalescenceModelType << endl << endl
            << "Valid coalescence models are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(logmom, dict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModel::coalescenceModel
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    logmom_(logmom)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModel::~coalescenceModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressible::momentumTransportModel&
Foam::coalescenceModel::continuousTurbulence() const
{
    return
        logmom_.phase().mesh()
       .lookupObject<phaseCompressible::momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                logmom_.continuousPhase().name()
            )
        );
}

// ************************************************************************* //

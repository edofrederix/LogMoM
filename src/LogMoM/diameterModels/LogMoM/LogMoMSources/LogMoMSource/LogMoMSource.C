#include "LogMoMSource.H"
#include "LogMoM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(LogMoMSource, 0);
    defineRunTimeSelectionTable(LogMoMSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::LogMoMSource>
Foam::diameterModels::LogMoMSource::New
(
    const word& type,
    const LogMoM& logmom,
    const dictionary& dict
)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown LogMoM source type "
            << type << nl << nl
            << "Valid LogMoM source types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<LogMoMSource>(cstrIter()(logmom, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ...

// ************************************************************************* //

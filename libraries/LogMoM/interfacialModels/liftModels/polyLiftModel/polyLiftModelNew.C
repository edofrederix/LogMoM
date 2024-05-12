#include "polyLiftModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyLiftModel> Foam::polyLiftModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer
)
{
    const dictionary& modelDict =
        outer ? interface.fluid().modelSubDict<polyLiftModel>(dict) : dict;

    const word polyLiftModelType(modelDict.lookup("type"));

    Info<< "Selecting polyLiftModel for "
        << interface.name() << ": " << polyLiftModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyLiftModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        return
            dictionaryConstructorTablePtr_->find("default")()
            (
                modelDict,
                interface
            );
    }
    else
    {
        return cstrIter()(modelDict, interface);
    }
}


Foam::autoPtr<Foam::blendedPolyLiftModel> Foam::blendedPolyLiftModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return autoPtr<blendedPolyLiftModel>(new blendedPolyLiftModel(dict, interface));
}


// ************************************************************************* //

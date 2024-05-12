#include "polyDragModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyDragModel> Foam::polyDragModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer,
    const bool registerObject
)
{
    const dictionary& modelDict =
        outer ? interface.fluid().modelSubDict<polyDragModel>(dict) : dict;

    const word polyDragModelType(modelDict.lookup("type"));

    Info<< "Selecting polyDragModel for "
        << interface.name() << ": " << polyDragModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyDragModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        return
            dictionaryConstructorTablePtr_->find("default")()
            (
                modelDict,
                interface,
                registerObject
            );
    }
    else
    {
        return cstrIter()(modelDict, interface, registerObject);
    }
}


Foam::autoPtr<Foam::blendedPolyDragModel> Foam::blendedPolyDragModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return autoPtr<blendedPolyDragModel>
    (
        new blendedPolyDragModel(dict, interface)
    );
}


// ************************************************************************* //

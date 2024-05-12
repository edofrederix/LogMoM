#include "polyTurbulentDispersionModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyTurbulentDispersionModel>
Foam::polyTurbulentDispersionModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer
)
{
    const dictionary& modelDict =
        outer
      ? interface.fluid().modelSubDict<polyTurbulentDispersionModel>(dict)
      : dict;

    const word polyTurbulentDispersionModelType(modelDict.lookup("type"));

    Info<< "Selecting polyTurbulentDispersionModel for "
        << interface.name() << ": " << polyTurbulentDispersionModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyTurbulentDispersionModelType);

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


Foam::autoPtr<Foam::blendedPolyTurbulentDispersionModel>
Foam::blendedPolyTurbulentDispersionModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return
        autoPtr<blendedPolyTurbulentDispersionModel>
        (
            new blendedPolyTurbulentDispersionModel(dict, interface)
        );
}



// ************************************************************************* //

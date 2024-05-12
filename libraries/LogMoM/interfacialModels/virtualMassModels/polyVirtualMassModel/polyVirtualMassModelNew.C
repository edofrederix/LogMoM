#include "polyVirtualMassModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyVirtualMassModel> Foam::polyVirtualMassModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer,
    const bool registerObject
)
{
    const dictionary& modelDict =
        outer
      ? interface.fluid().modelSubDict<polyVirtualMassModel>(dict)
      : dict;

    const word polyVirtualMassModelType(modelDict.lookup("type"));

    Info<< "Selecting polyVirtualMassModel for "
        << interface.name() << ": " << polyVirtualMassModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyVirtualMassModelType);

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


Foam::autoPtr<Foam::blendedPolyVirtualMassModel>
Foam::blendedPolyVirtualMassModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return
        autoPtr<blendedPolyVirtualMassModel>
        (
            new blendedPolyVirtualMassModel(dict, interface)
        );
}


// ************************************************************************* //

#include "polyWallLubricationModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyWallLubricationModel>
Foam::polyWallLubricationModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer
)
{
    const dictionary& modelDict =
        outer
      ? interface.fluid().modelSubDict<polyWallLubricationModel>(dict)
      : dict;

    const word polyWallLubricationModelType(modelDict.lookup("type"));

    Info<< "Selecting polyWallLubricationModel for "
        << interface.name() << ": " << polyWallLubricationModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyWallLubricationModelType);

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


Foam::autoPtr<Foam::blendedPolyWallLubricationModel>
Foam::blendedPolyWallLubricationModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return
        autoPtr<blendedPolyWallLubricationModel>
        (
            new blendedPolyWallLubricationModel(dict, interface)
        );
}


// ************************************************************************* //

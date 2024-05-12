#include "polyDragModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyDragModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(polyDragModel, 0);
    defineRunTimeSelectionTable(polyDragModel, dictionary);
}

const Foam::dimensionSet Foam::polyDragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyDragModel::polyDragModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, interface.name()),
            interface.mesh().time().timeName(),
            interface.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyDragModel::~polyDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyDragModel::writeData(Ostream& os) const
{
    return os.good();
}


Foam::tmp<Foam::volScalarField> Foam::blendedPolyDragModel::K
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyDragModel::K,
            "K",
            polyDragModel::dimK,
            false,
            gamma
        );
}


Foam::tmp<Foam::surfaceScalarField> Foam::blendedPolyDragModel::Kf
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyDragModel::Kf,
            "Kf",
            polyDragModel::dimK,
            false,
            gamma
        );
}


// ************************************************************************* //

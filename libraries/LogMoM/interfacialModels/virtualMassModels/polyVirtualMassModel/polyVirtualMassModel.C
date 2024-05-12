#include "polyVirtualMassModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyVirtualMassModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(polyVirtualMassModel, 0);
    defineRunTimeSelectionTable(polyVirtualMassModel, dictionary);
}

const Foam::dimensionSet Foam::polyVirtualMassModel::dimK(dimDensity);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyVirtualMassModel::polyVirtualMassModel
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

Foam::polyVirtualMassModel::~polyVirtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyVirtualMassModel::writeData(Ostream& os) const
{
    return os.good();
}


Foam::tmp<Foam::volScalarField> Foam::blendedPolyVirtualMassModel::K
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyVirtualMassModel::K,
            "K",
            polyVirtualMassModel::dimK,
            false,
            gamma
        );
}


Foam::tmp<Foam::surfaceScalarField> Foam::blendedPolyVirtualMassModel::Kf
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyVirtualMassModel::Kf,
            "Kf",
            polyVirtualMassModel::dimK,
            false,
            gamma
        );
}


// ************************************************************************* //

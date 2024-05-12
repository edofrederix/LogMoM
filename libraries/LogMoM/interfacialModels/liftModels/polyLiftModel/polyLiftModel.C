#include "polyLiftModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyLiftModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(polyLiftModel, 0);
    defineRunTimeSelectionTable(polyLiftModel, dictionary);
}

const Foam::dimensionSet Foam::polyLiftModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyLiftModel::polyLiftModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyLiftModel::~polyLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::blendedPolyLiftModel::F
(
    const label gamma
) const
{
    return evaluate(&polyLiftModel::F, "F", polyLiftModel::dimF, true, gamma);
}


Foam::tmp<Foam::surfaceScalarField> Foam::blendedPolyLiftModel::Ff
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyLiftModel::Ff,
            "Ff",
            polyLiftModel::dimF*dimArea,
            true,
            gamma
        );
}


// ************************************************************************* //

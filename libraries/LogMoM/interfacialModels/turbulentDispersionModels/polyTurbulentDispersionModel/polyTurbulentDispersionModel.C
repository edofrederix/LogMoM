#include "polyTurbulentDispersionModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyTurbulentDispersionModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug
    (
        polyTurbulentDispersionModel,
        0
    );
    defineRunTimeSelectionTable(polyTurbulentDispersionModel, dictionary);
}

const Foam::dimensionSet
Foam::polyTurbulentDispersionModel::dimD(1, -1, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyTurbulentDispersionModel::polyTurbulentDispersionModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyTurbulentDispersionModel::~polyTurbulentDispersionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendedPolyTurbulentDispersionModel::D
(
    const label gamma
) const
{
    return
        evaluate
        (
            &polyTurbulentDispersionModel::D,
            "F",
            polyTurbulentDispersionModel::dimD,
            true,
            gamma
        );
}


// ************************************************************************* //

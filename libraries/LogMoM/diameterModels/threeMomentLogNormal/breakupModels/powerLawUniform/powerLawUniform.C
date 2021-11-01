#include "powerLawUniform.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{
    defineTypeNameAndDebug(powerLawUniform, 0);
    addToRunTimeSelectionTable(breakupModel, powerLawUniform, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModels::powerLawUniform::powerLawUniform
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    breakupModel(pair, dict),
    power_
    (
        "power",
        dimless,
        dict
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::breakupModels::powerLawUniform::~powerLawUniform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::breakupModels::powerLawUniform::binaryRate
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    const scalar pi(constant::mathematical::pi);

    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            "tR",
            pair_.mesh(),
            dimensionedScalar(inv(dimTime*dimVolume), Zero)
        )
    );

    volScalarField& R = tR.ref();

    const volScalarField v1
    (
        pi*pow(d1,3)/6.0
    );

    const dimensioned<double> unitFix
    (
        dimensioned<double>(inv(pow(dimVolume, power_)*dimTime), 1.0)
    );

    R = 2.0*pow(v1, power_)*unitFix/v1;

    return tR;
}

// ************************************************************************* //

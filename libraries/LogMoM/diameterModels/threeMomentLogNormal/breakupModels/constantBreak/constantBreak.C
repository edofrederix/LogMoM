#include "constantBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{
    defineTypeNameAndDebug(constantBreak, 0);
    addToRunTimeSelectionTable(breakupModel, constantBreak, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModels::constantBreak::constantBreak
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    breakupModel(pair, dict),
    B_
    (
        "B",
        inv(dimTime),
        dict
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::breakupModels::constantBreak::~constantBreak()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::breakupModels::constantBreak::binaryRate
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    const scalar pi(constant::mathematical::pi);

    const volScalarField v1
    (
        pi*pow(d1,3)/6.0
    );

    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            "tR",
            pair_.mesh(),
            B_
        )
    );

    volScalarField& R = tR.ref();

    R *= 2.0/v1;

    return tR;
}

// ************************************************************************* //

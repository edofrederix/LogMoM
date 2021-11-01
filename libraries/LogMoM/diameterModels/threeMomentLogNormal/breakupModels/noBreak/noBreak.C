#include "noBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{
    defineTypeNameAndDebug(noBreak, 0);
    addToRunTimeSelectionTable(breakupModel, noBreak, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModels::noBreak::noBreak
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    breakupModel(pair, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::breakupModels::noBreak::~noBreak()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::breakupModels::noBreak::binaryRate
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            "tR",
            pair_.mesh(),
            dimensionedScalar(inv(dimTime), Zero)
        )
    );

    return tR;
}

// ************************************************************************* //

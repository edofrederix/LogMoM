#include "noBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

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
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    breakupModel(logmom, dict)
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
            logmom_.phase().mesh(),
            dimensionedScalar(inv(dimTime), Zero)
        )
    );

    return tR;
}

// ************************************************************************* //

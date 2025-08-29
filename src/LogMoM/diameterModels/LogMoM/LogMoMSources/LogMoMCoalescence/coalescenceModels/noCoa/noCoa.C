#include "noCoa.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(noCoa, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        noCoa,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::noCoa::noCoa
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    coalescenceModel(logmom, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModels::noCoa::~noCoa()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceModels::noCoa::rate
(
    const volScalarField& di,
    const volScalarField& dj
) const
{
    tmp<volScalarField> tE
    (
        volScalarField::New
        (
            "tE",
            logmom_.phase().mesh(),
            dimensionedScalar(dimVolume/dimTime, Zero)
        )
    );

    return tE;
}

// ************************************************************************* //

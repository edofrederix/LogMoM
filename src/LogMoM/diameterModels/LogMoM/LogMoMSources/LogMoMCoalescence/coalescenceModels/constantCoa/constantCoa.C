#include "constantCoa.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(constantCoa, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        constantCoa,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::constantCoa::constantCoa
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    coalescenceModel(logmom, dict),
    K_("K", dimVolume/dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModels::constantCoa::~constantCoa()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceModels::constantCoa::rate
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
            K_
        )
    );

    return tE;
}

// ************************************************************************* //

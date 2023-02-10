#include "noCoaFreq.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceFrequencyModels
{
    defineTypeNameAndDebug(noCoaFreq, 0);
    addToRunTimeSelectionTable(coalescenceFrequencyModel, noCoaFreq, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::noCoaFreq::noCoaFreq
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    coalescenceFrequencyModel(pair, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::noCoaFreq::~noCoaFreq()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyModels::noCoaFreq::frequency
(
    const volScalarField& di,
    const volScalarField& dj
) const
{
    tmp<volScalarField> tF
    (
        volScalarField::New
        (
            "tF",
            pair_.mesh(),
            dimensionedScalar(dimVolume/dimTime, Zero)
        )
    );

    return tF;
}

// ************************************************************************* //

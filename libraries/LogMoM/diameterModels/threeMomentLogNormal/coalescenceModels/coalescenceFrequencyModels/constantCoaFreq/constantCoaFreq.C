#include "constantCoaFreq.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceFrequencyModels
{
    defineTypeNameAndDebug(constantCoaFreq, 0);
    addToRunTimeSelectionTable(coalescenceFrequencyModel, constantCoaFreq, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::constantCoaFreq::constantCoaFreq
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    coalescenceFrequencyModel(pair, dict),
    K_
    (
        "K",
        dimVolume/dimTime,
        dict.subDict("constantFrequencyCoeffs")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::constantCoaFreq::~constantCoaFreq()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyModels::constantCoaFreq::frequency
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
            K_
        )
    );

    return tF;
}

// ************************************************************************* //

#include "constantCoaFreq.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

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
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    coalescenceFrequencyModel
    (
        pair,
        dict.subDict(this->type() + this->coeffsDictName_)
    ),
    K_("K", dimVolume/dimTime, coeffs_)
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

#include "LehrMilliesMewesCoaFreq.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "phaseSystem.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceFrequencyModels
{
    defineTypeNameAndDebug(LehrMilliesMewesCoaFreq, 0);
    addToRunTimeSelectionTable
    (
        coalescenceFrequencyModel,
        LehrMilliesMewesCoaFreq,
        dictionary
    );
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::LehrMilliesMewesCoaFreq
::LehrMilliesMewesCoaFreq
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    coalescenceFrequencyModel(pair, dict),
    uCrit_
    (
        dimensionedScalar::lookupOrDefault("uCrit", dict, dimVelocity, 0.08)
    ),
    alphaMax_
    (
        dimensionedScalar::lookupOrDefault("alphaMax", dict, dimless, 0.6)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::LehrMilliesMewesCoaFreq
::~LehrMilliesMewesCoaFreq()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyModels::LehrMilliesMewesCoaFreq::frequency
(
    const volScalarField& d1,
    const volScalarField& d2
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

    volScalarField& F = tF.ref();

    const scalar pi(constant::mathematical::pi);
    const phaseModel& contPhase = pair_.continuous();


    const volScalarField uChar
    (
        sqrt(2.0)*cbrt(continuousTurbulence().epsilon())
       *sqrt(cbrt(sqr(d1)) + cbrt(sqr(d2)))
    );

    // This formulation includes the coalescence efficiency. Therefore, a
    // constant efficiency coalescence with K=1 should be used.

    F +=
      pi/4.0*sqr(d1 + d2)*min(uChar, uCrit_)
     *exp
      (
        - sqr(cbrt(alphaMax_)
         /cbrt(max(contPhase, 0.0001)) - 1.0)
      );


    return tF;
}

// ************************************************************************* //

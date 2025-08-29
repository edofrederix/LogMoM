#include "LehrMilliesMewesCoa.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "phaseSystem.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(LehrMilliesMewes, 0);
    addToRunTimeSelectionTable(coalescenceModel,LehrMilliesMewes,dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::LehrMilliesMewes::LehrMilliesMewes
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    coalescenceModel(logmom, dict),
    uCrit_("uCrit", dimVelocity, dict, 0.08),
    alphaMax_("alphaMax", dimless, dict, 0.6)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModels::LehrMilliesMewes::~LehrMilliesMewes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceModels::LehrMilliesMewes::rate
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    const scalar pi(constant::mathematical::pi);

    tmp<volScalarField> tF
    (
        volScalarField::New
        (
            "tF",
            logmom_.phase().mesh(),
            dimensionedScalar(dimVolume/dimTime, Zero)
        )
    );

    volScalarField& F = tF.ref();

    const phaseModel& phase = logmom_.interface().dispersed();
    const volScalarField epsilon(continuousTurbulence().epsilon());

    const volScalarField uPrime
    (
        max
        (
            sqrt(2.0)*cbrt(epsilon)*sqrt(cbrt(sqr(d1)) + cbrt(sqr(d2))),
            dimensionedScalar(dimVelocity, 0)
        )
    );

    // Eq. (20) of Lehr et al. (2002)

    F =
        pi/4.0
      * sqr(d1 + d2)
      * min(uPrime, uCrit_)
      * exp
        (
          - sqr
            (
                cbrt(alphaMax_/max(phase, phase.residualAlpha()))
              - 1.0
            )
        );

    return tF;
}

// ************************************************************************* //

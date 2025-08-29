#include "LehrMilliesMewesBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "phaseSystem.H"
#include "mathematicalConstants.H"
#include "linearInterpolationWeights.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{
    defineTypeNameAndDebug(LehrMilliesMewes, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        LehrMilliesMewes,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModels::LehrMilliesMewes::LehrMilliesMewes
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    breakupModel(logmom, dict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::breakupModels::LehrMilliesMewes::binaryRate
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    const scalar pi(constant::mathematical::pi);

    tmp<volScalarField> tR
    (
        volScalarField::New
        (
            "tR",
            logmom_.phase().mesh(),
            dimensionedScalar(inv(dimTime*dimVolume), Zero)
        )
    );

    volScalarField& R = tR.ref();

    const phaseModel& continuousPhase = logmom_.continuousPhase();
    const volScalarField epsilon(continuousTurbulence().epsilon());

    const dimensionedScalar epsilonMin(epsilon.dimensions(), SMALL);

    const volScalarField L
    (
        pow(logmom_.interface().sigma()/continuousPhase.rho(), 3.0/5.0)
      / pow(max(epsilon,epsilonMin), 2.0/5.0)
    );

    const volScalarField T
    (
        pow(logmom_.interface().sigma()/continuousPhase.rho(), 2.0/5.0)
      / pow(max(epsilon,epsilonMin), 3.0/5.0)
    );

    // When the volume of the daughter is larger than half the parent volume,
    // use symmetry about half the parent volume (see Eq. (16) and (17) of Lehr
    // et al. (2002))

    const volScalarField f(neg0(d2 - cbrt(0.5)*d1));
    const volScalarField d2p
    (
        f*d2 + (1.0-f)*cbrt(pow3(d1) - pow3(d2))
    );

    // Eq. (12-17) of Lehr et al. (2002)

    R =
        0.5*pow(d1/L, 5.0/3.0)
      * exp(-sqrt(2.0)/pow3(d1/L))
      * 6.0/pow(pi, 1.5)/pow3(d2p/L)
      * exp(-9.0/4.0*sqr(log(pow(2.0, 0.4)*d2p/L)))
      / max(1.0 + erf(1.5*log(pow(2.0, 1.0/15.0)*d1/L)), SMALL)
      / (T*pow3(L));

    return tR;
}

// ************************************************************************* //

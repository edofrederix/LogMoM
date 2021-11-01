#include "LehrMilliesMewesBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "phaseSystem.H"
#include "mathematicalConstants.H"
#include "linearInterpolationWeights.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{
    defineTypeNameAndDebug(LehrMilliesMewesBreak, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        LehrMilliesMewesBreak,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModels::LehrMilliesMewesBreak::LehrMilliesMewesBreak
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    breakupModel(pair, dict),
    sigma_
    (
        "sigma",
        dimMass/sqr(dimTime),
        dict
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::breakupModels::LehrMilliesMewesBreak::binaryRate
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
            dimensionedScalar(inv(dimTime*dimVolume), Zero)
        )
    );

    volScalarField& R = tR.ref();

    const phaseModel& continuousPhase = pair_.continuous();
    const volScalarField epsilon(continuousTurbulence().epsilon());

    volScalarField L
    (
        pow
        (
            sigma_/continuousPhase.rho(),
            3.0/5.0
        )
       /pow(epsilon, 2.0/5.0)
    );

    // Reset of dimension to pure length to avoid problems in transcendental
    // functions due to small exponents

    L.dimensions().reset(dimLength);

    const volScalarField T
    (
        pow
        (
            sigma_/continuousPhase.rho(),
            2.0/5.0
        )
       /pow(epsilon, 3.0/5.0)
    );

    const scalar pi(constant::mathematical::pi);

    if (d2<pow(0.5, 1.0/3.0)*d1)
    {
        R = 0.5*pow(d1/L, 5.0/3.0)
          * exp(-sqrt(2.0)/pow3(d1/L))
          * 6.0/pow(pi, 1.5)
          / pow3
            (
                max
                (
                    min(d2, pow(pow3(d1)*1.0, 1.0/3.0)),
                    pow(pow3(d1)*0.0, 1.0/3.0)
                )
              / L
            )
          * exp
            (
              - 9.0/4.0
              * sqr
                (
                    log
                    (
                        pow(2.0, 0.4)
                      * max
                          (
                              min(d2, pow(pow3(d1)*1.0, 1.0/3.0)),
                              pow(pow3(d1)*0.0, 1.0/3.0)
                          )
                      / L
                    )
                )
            )
          / max(1.0 + erf(1.5*log(pow(2.0, 1.0/15.0)*d1/L)), SMALL)
          / (T*pow3(L));
    }
    else
    {
        R =  0.5*pow(d1/L, 5.0/3.0)
          * exp(-sqrt(2.0)/pow3(d1/L))
          * 6.0/pow(pi, 1.5)
          / pow3
            (
                pow
                (
                    pow3(d1)
                  - pow3
                    (
                        max
                        (
                            min(d2, pow(pow3(d1)*1.0, 1.0/3.0)),
                            pow(pow3(d1)*0.0, 1.0/3.0)
                        )
                    ),
                    1.0/3.0
                )
              / L
            )
          * exp
            (
              - 9.0/4.0
              * sqr
                (
                    log
                    (
                        pow(2.0, 0.4)
                      * pow
                        (
                            pow3(d1)
                          - pow3
                            (
                                max
                                (
                                    min(d2, pow(pow3(d1)*1.0, 1.0/3.0)),
                                    pow(pow3(d1)*0.0, 1.0/3.0)
                                )
                            ),
                            1.0/3.0
                        )
                      / L
                    )
                )
            )
          / max(1.0 + erf(1.5*log(pow(2.0, 1.0/15.0)*d1/L)), SMALL)
          / (T*pow3(L));
    }

    return tR;
}

// ************************************************************************* //

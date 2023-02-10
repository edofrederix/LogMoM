#include "LuoSvendsenBreak.H"
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
    defineTypeNameAndDebug(LuoSvendsenBreak, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        LuoSvendsenBreak,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModels::LuoSvendsenBreak::LuoSvendsenBreak
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    breakupModel(pair, dict),
    gammaUpperReg5by11_(),
    C4_(dimensionedScalar::lookupOrDefault("C4", dict, dimless, 0.923)),
    beta_(dimensionedScalar::lookupOrDefault("beta", dict, dimless, 2.05)),
    minEddyRatio_
    (
        dimensionedScalar::lookupOrDefault
        (
            "minEddyRatio",
            dict,
            dimless,
            11.4
        )
    ),
    sigma_
    (
        "sigma",
        dimMass/sqr(dimTime),
        dict
    )
{
    List<Tuple2<scalar, scalar>> gammaUpperReg5by11Table;

    gammaUpperReg5by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));

    for (scalar z = 1e-2; z <= 10.0; z = z + 1e-2)
    {
        Tuple2<scalar, scalar> gamma5by11
        (
            z,
            incGammaRatio_Q(5.0/11.0, z)
        );

        gammaUpperReg5by11Table.append(gamma5by11);
    }

    gammaUpperReg5by11_ =
        new Function1s::Table<scalar>
        (
            "gamma5by11",
            Function1s::tableBase::boundsHandling::clamp,
            linearInterpolationWeights::typeName,
            autoPtr<TableReader<scalar>>(nullptr),
            gammaUpperReg5by11Table
        );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::breakupModels::LuoSvendsenBreak::binaryRate
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

    const scalar pi(constant::mathematical::pi);

    const volScalarField epsilon(continuousTurbulence().epsilon());

    const phaseModel& contPhase = pair_.continuous();

    const volScalarField eta
    (
        pow025(pow3(contPhase.thermo().nu())/epsilon)
    );

    const volScalarField rho(contPhase.rho());

    const volScalarField xiMin(minEddyRatio_*eta/d1);

    const volScalarField v1(pi*pow(d1, 3.0)/6.0);
    const volScalarField v2(pi*pow(d2, 3.0)/6.0);

    const volScalarField fBV(max(min(v2/v1, 0.95), 0.05));

    const volScalarField cf
    (
        pow(fBV,2.0/3.0) + pow(1.0-fBV,2.0/3.0) - 1.0
    );

    const volScalarField b
    (
        12.0*cf*sigma_
      / (
            beta_*rho
          * pow(epsilon, 2.0/3.0)
          * pow(d1, 5.0/3.0)
        )
    );

    const volScalarField bMin(b/pow(xiMin, 11.0/3.0));

    volScalarField integral(3.0/(11.0*pow(b,8.0/11.0)));

    forAll(integral, celli)
    {
        integral[celli] *=
            2.0*pow(b[celli], 3.0/11.0)*tgamma(5.0/11.0)
          * (
                gammaUpperReg5by11_->value(b[celli])
              - gammaUpperReg5by11_->value(bMin[celli])
            );
    }

    R =
        max
        (
            C4_*contPhase/v1
          * cbrt(epsilon/sqr(d1))
          * integral,
            dimensionedScalar
            (
                inv(dimTime*dimVolume),
                Zero
            )
        );

    return tR;
}

// ************************************************************************* //

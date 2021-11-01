#include "PrinceBlanchCoaFreq.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "phaseSystem.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceFrequencyModels
{
    defineTypeNameAndDebug(PrinceBlanchCoaFreq, 0);
    addToRunTimeSelectionTable
    (
        coalescenceFrequencyModel,
        PrinceBlanchCoaFreq,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::PrinceBlanchCoaFreq::PrinceBlanchCoaFreq
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    coalescenceFrequencyModel(pair, dict),
    C1_(dimensionedScalar::lookupOrDefault("C1", dict, dimless, 0.356)),
    sigma_("sigma", dimMass/sqr(dimTime), dict.lookup("sigma")),
    turbulent_(dict.lookup("turbulentCoalescence")),
    buoyant_(dict.lookup("buoyantCoalescence")),
    laminar_(dict.lookup("laminarCoalescence"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::PrinceBlanchCoaFreq::~PrinceBlanchCoaFreq()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyModels::PrinceBlanchCoaFreq::frequency
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

    if (turbulent_)
    {
        const phaseCompressible::momentumTransportModel& turbulence =
            continuousTurbulence();

        F +=
            C1_*pi*sqr(d1+d2)
          * cbrt(turbulence.epsilon())
          * sqrt(pow(d1,2.0/3.0)+pow(d2,2.0/3.0));
    }

    if (buoyant_)
    {
        const uniformDimensionedVectorField& g =
            pair_.mesh().lookupObject<uniformDimensionedVectorField>("g");

        const volScalarField Sij(pi/4.0*sqr(d1+d2));

        F +=
            Sij
          * mag
            (
                sqrt
                (
                    2.14*sigma_
                  / (pair_.continuous().rho()*d1)
                  + 0.505*mag(g)*d1
                )
              - sqrt
                (
                    2.14*sigma_
                  / (pair_.continuous().rho()*d2)
                  + 0.505*mag(g)*d2
                )
            );

    }

    if (laminar_)
    {
        const volScalarField shearStrainRate
        (
            sqrt(2.0)*mag(symm(fvc::grad(pair_.continuous().U())))
        );

        F += 1.0/6.0*pow3(d1+d2)*shearStrainRate;
    }

    return tF;
}

// ************************************************************************* //

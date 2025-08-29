#include "PrinceBlanch.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "phaseSystem.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(PrinceBlanch, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        PrinceBlanch,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::PrinceBlanch::PrinceBlanch
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    coalescenceModel(logmom, dict),
    h0_("h0", dimLength, dict, 1e-4),
    hf_("hf", dimLength, dict, 1e-8),
    C1_("C1", dimless, dict, 0.356),
    turbulence_(dict.lookup("turbulence")),
    buoyancy_(dict.lookup("buoyancy")),
    laminarShear_(dict.lookup("laminarShear"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceModels::PrinceBlanch::~PrinceBlanch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceModels::PrinceBlanch::rate
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    const scalar pi(constant::mathematical::pi);

    // Efficiency

    tmp<volScalarField> tE
    (
        volScalarField::New
        (
            "E",
            logmom_.phase().mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );

    volScalarField& E = tE.ref();

    const phaseCompressible::momentumTransportModel& turbulence =
        continuousTurbulence();

    const volScalarField rij(1.0/(1.0/d1 + 1.0/d2));

    E =
        exp
        (
          - sqrt
            (
                pow3(rij)*logmom_.continuousPhase().rho()
              / (16.0*logmom_.interface().sigma())
            )
          * log(h0_/hf_)
          * cbrt(turbulence.epsilon())
          / pow(rij,2.0/3.0)
        );

    // Frequency

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

    if (turbulence_)
    {
        F +=
            C1_*pi*sqr(d1+d2)
          * cbrt(turbulence.epsilon())
          * sqrt(pow(d1,2.0/3.0)+pow(d2,2.0/3.0));
    }

    if (buoyancy_)
    {
        const uniformDimensionedVectorField& g =
            logmom_.phase().mesh()
           .lookupObject<uniformDimensionedVectorField>("g");

        const volScalarField Sij(pi/4.0*sqr(d1+d2));

        F +=
            Sij
          * mag
            (
                sqrt
                (
                    2.14*logmom_.interface().sigma()
                  / (logmom_.continuousPhase().rho()*d1)
                  + 0.505*mag(g)*d1
                )
              - sqrt
                (
                    2.14*logmom_.interface().sigma()
                  / (logmom_.continuousPhase().rho()*d2)
                  + 0.505*mag(g)*d2
                )
            );

    }

    if (laminarShear_)
    {
        const volScalarField shearStrainRate
        (
            sqrt(2.0)*mag(symm(fvc::grad(logmom_.continuousPhase().U())))
        );

        F += 1.0/6.0*pow3(d1+d2)*shearStrainRate;
    }

    return E*F;
}

// ************************************************************************* //

#include "PrinceBlanchCoaEff.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceEfficiencyModels
{
    defineTypeNameAndDebug(PrinceBlanchCoaEff, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyModel,
        PrinceBlanchCoaEff,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::PrinceBlanchCoaEff::PrinceBlanchCoaEff
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    coalescenceEfficiencyModel(pair, dict),
    h0_
    (
        dimensionedScalar::lookupOrDefault
        (
            "h0",
            dict,
            dimLength,
            1e-4
        )
    ),
    hf_
    (
        dimensionedScalar::lookupOrDefault
        (
            "hf",
            dict,
            dimLength,
            1e-8
        )
    ),
    sigma_("sigma", dimMass/sqr(dimTime), dict.lookup("sigma"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::PrinceBlanchCoaEff::~PrinceBlanchCoaEff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyModels::PrinceBlanchCoaEff::efficiency
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{

    tmp<volScalarField> tE
    (
        volScalarField::New
        (
            "tE",
            pair_.mesh(),
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
                pow3(rij)*pair_.continuous().rho()
              / (16.0*sigma_)
            )
          * log(h0_/hf_)
          * cbrt(turbulence.epsilon())
          / pow(rij,2.0/3.0)
        );

    return tE;
}

// ************************************************************************* //

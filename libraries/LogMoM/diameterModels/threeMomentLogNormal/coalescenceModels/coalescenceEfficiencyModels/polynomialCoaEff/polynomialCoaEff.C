#include "polynomialCoaEff.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceEfficiencyModels
{
    defineTypeNameAndDebug(polynomialCoaEff, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyModel,
        polynomialCoaEff,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::polynomialCoaEff::polynomialCoaEff
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    coalescenceEfficiencyModel(pair, dict),
    K_(dict.subDict("polynomialEfficiencyCoeffs").lookup("K")),
    p_(dict.subDict("polynomialEfficiencyCoeffs").lookup("p")),
    q_(dict.subDict("polynomialEfficiencyCoeffs").lookup("q"))
{
    if (K_.size() != p_.size() || K_.size() != q_.size())
    {
        FatalErrorInFunction
            << "Lists K, p and q must have equal length"
            << abort(FatalError);
    }

    if (K_.size() == 0)
    {
        FatalErrorInFunction
            << "Lists K, p and q should have at least one entry"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::polynomialCoaEff::~polynomialCoaEff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyModels::polynomialCoaEff::efficiency
(
    const volScalarField& di,
    const volScalarField& dj
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

    for (label i = 0; i < K_.size(); i++)
    {
        const scalar p = p_[i];
        const scalar q = q_[i];

        const dimensionedScalar K(E.dimensions()/pow(dimLength,p+q), K_[i]);

        E += K*(pow(di,p)*pow(dj,q) + pow(di,q)*pow(dj,p));
    }

    return tE;
}

// ************************************************************************* //

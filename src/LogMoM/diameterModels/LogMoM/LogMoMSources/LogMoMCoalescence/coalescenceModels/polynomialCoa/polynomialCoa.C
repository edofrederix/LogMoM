#include "polynomialCoa.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(polynomialCoa, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        polynomialCoa,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceModels::polynomialCoa::polynomialCoa
(
    const diameterModels::LogMoM& logmom,
    const dictionary& dict
)
:
    coalescenceModel(logmom, dict),
    K_(dict.lookup("K")),
    p_(dict.lookup("p")),
    q_(dict.lookup("q"))
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

Foam::coalescenceModels::polynomialCoa::~polynomialCoa()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceModels::polynomialCoa::rate
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
            logmom_.phase().mesh(),
            dimensionedScalar(dimVolume/dimTime, Zero)
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

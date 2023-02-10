#include "polynomialCoaFreq.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceFrequencyModels
{
    defineTypeNameAndDebug(polynomialCoaFreq, 0);
    addToRunTimeSelectionTable
    (
        coalescenceFrequencyModel,
        polynomialCoaFreq,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyModels::polynomialCoaFreq::polynomialCoaFreq
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    coalescenceFrequencyModel(pair, dict),
    K_(dict.subDict("polynomialFrequencyCoeffs").lookup("K")),
    p_(dict.subDict("polynomialFrequencyCoeffs").lookup("p")),
    q_(dict.subDict("polynomialFrequencyCoeffs").lookup("q"))
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

Foam::coalescenceFrequencyModels::polynomialCoaFreq::~polynomialCoaFreq()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyModels::polynomialCoaFreq::frequency
(
    const volScalarField& di,
    const volScalarField& dj
) const
{
    tmp<volScalarField> tF
    (
        volScalarField::New
        (
            "tF",
            pair_.mesh(),
            dimVolume/dimTime
        )
    );

    volScalarField& F = tF.ref();

    for (label i = 0; i < K_.size(); i++)
    {
        const scalar p(p_[i]);
        const scalar q(q_[i]);

        const dimensionedScalar K(F.dimensions()/pow(dimLength,p+q), K_[i]);

        F += K*(pow(di,p)*pow(dj,q) + pow(di,q)*pow(dj,p));
    }

    return tF;
}

// ************************************************************************* //

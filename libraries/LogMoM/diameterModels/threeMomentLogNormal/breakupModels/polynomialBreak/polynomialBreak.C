#include "polynomialBreak.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{
    defineTypeNameAndDebug(polynomialBreak, 0);
    addToRunTimeSelectionTable(breakupModel, polynomialBreak, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::breakupModels::polynomialBreak::polynomialBreak
(
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    breakupModel
    (
        pair,
        dict.subDict(this->type() + this->coeffsDictName_)
    ),
    B_(coeffs_.lookup("B")),
    p_(coeffs_.lookup("p"))
{
    if (B_.size() != p_.size())
    {
        FatalErrorInFunction
            << "Lists B and p must have equal length"
            << abort(FatalError);
    }

    if (B_.size() == 0)
    {
        FatalErrorInFunction
            << "Lists B and p should have at least one entry"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::breakupModels::polynomialBreak::~polynomialBreak()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::breakupModels::polynomialBreak::binaryRate
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    const scalar pi(constant::mathematical::pi);

    const volScalarField v1
    (
        pi*pow(d1,3)/6.0
    );

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

    for (label i = 0; i < B_.size(); i++)
    {
        const scalar p(p_[i]);

        const dimensionedScalar B
        (
            R.dimensions()/pow(dimLength,p)*dimVolume,
            B_[i]
        );

        R += 2.0/v1*B*pow(d1,p);
    }

    return tR;
}

// ************************************************************************* //

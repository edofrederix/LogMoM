#include "inletOutletLogNormalFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::inletOutletLogNormalFvPatchScalarField::computeField() const
{
    const scalar pi(constant::mathematical::pi);

    const dimensionSet dims(internalField().dimensions());

    if
    (
        dims[dimensionSet::MASS] != 0.0
     || dims[dimensionSet::TIME] != 0.0
     || dims[dimensionSet::TEMPERATURE] != 0.0
     || dims[dimensionSet::MOLES] != 0.0
     || dims[dimensionSet::CURRENT] != 0.0
     || dims[dimensionSet::LUMINOUS_INTENSITY] != 0.0
    )
    {
        FatalErrorInFunction
            << "Dimensions of this moment field are incorrect"
            << abort(FatalError);
    }

    const scalar gamma(dims[dimensionSet::LENGTH]+3);

    scalar q(1.0);

    if (gamma == 0.0)
    {
        q = 1e-6;
    }
    else if (gamma == 2.0)
    {
        q = pi;
    }
    else if (gamma == 3.0)
    {
        q = pi/6.0;
    }

    return
        tmp<scalarField>
        (
            new scalarField
            (
                patch().size(),
                q*6.0/pi
              * pow(dsm_, gamma-3.0)
              * exp
                (
                    (0.5*sqr(gamma)-2.5*gamma+3.0)*sqr(sigma_)
                )
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletOutletLogNormalFvPatchScalarField::
inletOutletLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    sigma_(0.0),
    dsm_(0.0),
    alphaName_("")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::inletOutletLogNormalFvPatchScalarField::
inletOutletLogNormalFvPatchScalarField
(
    const inletOutletLogNormalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    sigma_(ptf.sigma_),
    dsm_(ptf.dsm_),
    alphaName_(ptf.alphaName_)
{}


Foam::inletOutletLogNormalFvPatchScalarField::
inletOutletLogNormalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    sigma_(readScalar(dict.lookup("sigma"))),
    dsm_(readScalar(dict.lookup("dsm"))),
    alphaName_(dict.lookup("alphaName"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    this->refValue() = Zero;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(0.0);
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}



Foam::inletOutletLogNormalFvPatchScalarField::
inletOutletLogNormalFvPatchScalarField
(
    const inletOutletLogNormalFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    sigma_(tppsf.sigma_),
    dsm_(tppsf.dsm_),
    alphaName_(tppsf.alphaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inletOutletLogNormalFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);
}


void Foam::inletOutletLogNormalFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(ptf, addr);
}


void Foam::inletOutletLogNormalFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = computeField();
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::inletOutletLogNormalFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);

    writeEntry(os, "sigma", sigma_);
    writeEntry(os, "dsm", dsm_);
    writeEntry(os, "alphaName", alphaName_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inletOutletLogNormalFvPatchScalarField
    );
}

// ************************************************************************* //

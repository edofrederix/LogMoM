#include "oscillatingDiskFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oscillatingDiskFvPatchScalarField::
oscillatingDiskFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    lambda_(dict.lookup<scalar>("lambda", dimLength)),
    kappa_(dict.lookup<scalar>("kappa", dimTime)),
    average_(dict.lookup<scalar>("average", iF.dimensions())),
    amplitude_(dict.lookup<scalar>("amplitude", iF.dimensions()))
{
    fvPatchScalarField::operator==(waves(this->db().time().value()));
}


Foam::oscillatingDiskFvPatchScalarField::
oscillatingDiskFvPatchScalarField
(
    const oscillatingDiskFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper, false), // Don't map
    lambda_(ptf.lambda_),
    kappa_(ptf.kappa_),
    average_(ptf.average_),
    amplitude_(ptf.amplitude_)
{
    // Set the patch pressure to the current total pressure
    // This is not ideal but avoids problems with the creation of patch faces
    fvPatchScalarField::operator==(waves(this->db().time().value()));
}


Foam::oscillatingDiskFvPatchScalarField::
oscillatingDiskFvPatchScalarField
(
    const oscillatingDiskFvPatchScalarField& ptf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    lambda_(ptf.lambda_),
    kappa_(ptf.kappa_),
    average_(ptf.average_),
    amplitude_(ptf.amplitude_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::oscillatingDiskFvPatchScalarField::waves(const scalar t) const
{
    const scalar pi(constant::mathematical::pi);

    const vectorField& Cf = patch().Cf();
    const scalarField& magSf = patch().magSf();

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const vector center = gSum(magSf*Cf)/gSum(magSf);

    tmp<scalarField> tValues(new scalarField(patch().size()));
    scalarField& values = tValues.ref();

    const scalarField r(Foam::mag(Cf - center));
    values = average_ + amplitude_*cos(2.0*pi/lambda_*r)*cos(2.0*pi/kappa_*t);

    // Restore tag
    UPstream::msgType() = oldTag;

    return tValues;
}

void Foam::oscillatingDiskFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchScalarField::operator==(waves(this->db().time().value()));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::oscillatingDiskFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "lambda", lambda_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "average", average_);
    writeEntry(os, "amplitude", amplitude_);

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        oscillatingDiskFvPatchScalarField
    );
}

// ************************************************************************* //

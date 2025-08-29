#include "TOPFLOWAlphaInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TOPFLOWAlphaInletFvPatchScalarField::
TOPFLOWAlphaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    average_(0.0),
    phaseName_("air"),
    caseName_("A")
{}


Foam::TOPFLOWAlphaInletFvPatchScalarField::
TOPFLOWAlphaInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    average_(dict.lookup<scalar>("average")),
    phaseName_(dict.lookup<word>("phase")),
    caseName_(dict.lookup<word>("case"))
{
    if (caseName_ != "A" && caseName_ != "B" && caseName_ != "C")
    {
        FatalErrorInFunction
            << "Invalid TOPFLOW case specified"
            << endl << abort(FatalError);
    }

    if (phaseName_ != "air" && phaseName_ != "water")
    {
        FatalErrorInFunction
            << "Invalid phase specified"
            << endl << abort(FatalError);
    }

    IFstream file("data/TOPFLOW_" + caseName_ + ".txt");
    file >> data_;

    fixedValueFvPatchScalarField::evaluate();
}


Foam::TOPFLOWAlphaInletFvPatchScalarField::
TOPFLOWAlphaInletFvPatchScalarField
(
    const TOPFLOWAlphaInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    average_(ptf.average_),
    phaseName_(ptf.phaseName_),
    caseName_(ptf.caseName_),
    data_(ptf.data_)
{}


Foam::TOPFLOWAlphaInletFvPatchScalarField::
TOPFLOWAlphaInletFvPatchScalarField
(
    const TOPFLOWAlphaInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    average_(ptf.average_),
    phaseName_(ptf.phaseName_),
    caseName_(ptf.caseName_),
    data_(ptf.data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TOPFLOWAlphaInletFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarList r(data_.size());
    scalarList alpha(data_.size());

    const vectorField& Cf = patch().Cf();
    const scalarField& magSf = patch().magSf();

    forAll(data_, i)
    {
        r[i] = data_[i].first();
        alpha[i] = max(data_[i].second(),1e-5);
    }

    scalarField alphaPatchAir(Cf.size());

    forAll(Cf, facei)
    {
        const scalar ri = Cf[facei].x();

        const label i = findLower(r,ri);

        if (i == -1)
        {
            alphaPatchAir[facei] = alpha[0];
        }
        else if (i < r.size()-1)
        {
            alphaPatchAir[facei] =
                (ri     - r[i])/(r[i+1] - r[i]) * alpha[i]
              + (r[i+1] - ri  )/(r[i+1] - r[i]) * alpha[i+1];
        }
        else
        {
            alphaPatchAir[facei] = alpha[alpha.size()-1];
        }
    }

    // Normalize

    alphaPatchAir =
        alphaPatchAir*gSum(magSf)/gSum(alphaPatchAir*magSf);

    // Set

    if (phaseName_ == "air")
    {
        fixedValueFvPatchScalarField::operator==
        (
            alphaPatchAir*average_
        );
    }
    else
    {
        fixedValueFvPatchScalarField::operator==
        (
            1.0 - alphaPatchAir*(1.0-average_)
        );
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::TOPFLOWAlphaInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "average", average_);
    writeEntry(os, "phase", phaseName_);
    writeEntry(os, "case", caseName_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        TOPFLOWAlphaInletFvPatchScalarField
    );
}

// ************************************************************************* //

#include "inletVelocityMappedInternalValueFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cell_interpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletVelocityMappedInternalValueFvPatchVectorField::
inletVelocityMappedInternalValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    mappedInternalValueFvPatchField<vector>(p, iF, dict)
{}

Foam::inletVelocityMappedInternalValueFvPatchVectorField::
inletVelocityMappedInternalValueFvPatchVectorField
(
    const inletVelocityMappedInternalValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    mappedInternalValueFvPatchField<vector>(ptf, p, iF, mapper)
{}

Foam::inletVelocityMappedInternalValueFvPatchVectorField::
inletVelocityMappedInternalValueFvPatchVectorField
(
    const inletVelocityMappedInternalValueFvPatchVectorField& ptf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    mappedInternalValueFvPatchField<vector>(ptf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inletVelocityMappedInternalValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const fvMesh& nbrMesh = refCast<const fvMesh>(this->mapper().nbrMesh());

    const VolField<vector>& nbrField =
        this->mapper().sameRegion()
     && this->fieldName_ == this->internalField().name()
      ? refCast<const VolField<vector>>(this->internalField())
      : nbrMesh.template lookupObject<VolField<vector>>(this->fieldName_);

    // Construct mapped values
    Field<vector> sampleValues;

    if (this->interpolationScheme_ != interpolations::cell<vector>::typeName)
    {
        // Create an interpolation
        autoPtr<interpolation<vector>> interpolatorPtr
        (
            interpolation<vector>::New
            (
                this->interpolationScheme_,
                nbrField
            )
        );
        const interpolation<vector>& interpolator = interpolatorPtr();

        // Cells on which samples are generated
        const labelList& sampleCells = this->mapper().cellIndices();

        // Send the patch points to the cells
        pointField samplePoints(this->mapper().samplePoints());
        this->mapper().map().reverseDistribute
        (
            sampleCells.size(),
            samplePoints
        );

        // Interpolate values
        sampleValues.resize(sampleCells.size());
        forAll(sampleCells, i)
        {
            if (sampleCells[i] != -1)
            {
                sampleValues[i] =
                    interpolator.interpolate
                    (
                        samplePoints[i],
                        sampleCells[i]
                    );
            }
        }

        // Send the values back to the patch
        this->mapper().map().distribute(sampleValues);
    }
    else
    {
        // No interpolation. Just sample cell values directly.
        sampleValues = this->mapper().distribute(nbrField);
    }

    // Strip tangential components and limit outflow

    sampleValues =
        Foam::min(this->patch().nf() & sampleValues, 0.0)*this->patch().nf();

    // Set the average, if necessary
    if (this->setAverage_)
    {

        const vector sampleAverageValue =
            gSum(this->patch().magSf()*sampleValues)
           /gSum(this->patch().magSf());

        if (mag(sampleAverageValue)/mag(this->average_) > 0.5)
        {
            sampleValues *= mag(this->average_)/mag(sampleAverageValue);
        }
        else
        {
            sampleValues += this->average_ - sampleAverageValue;
        }
    }

    // Assign sampled patch values
    this->operator==(sampleValues);

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchField<vector>::updateCoeffs();
}

void Foam::inletVelocityMappedInternalValueFvPatchVectorField::write
(
    Ostream& os
) const
{
    mappedInternalValueFvPatchField<vector>::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        inletVelocityMappedInternalValueFvPatchVectorField
    );
}

// ************************************************************************* //

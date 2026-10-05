#include "weightedMappedInternalValueFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cell_interpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::weightedMappedInternalValueFvPatchField<Type>::
weightedMappedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, fvMesh>& iF,
    const dictionary& dict
)
:
    mappedInternalValueFvPatchField<Type>(p, iF, dict),
    phiName_(dict.lookup<word>("phi"))
{}

template<class Type>
Foam::weightedMappedInternalValueFvPatchField<Type>::
weightedMappedInternalValueFvPatchField
(
    const weightedMappedInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    mappedInternalValueFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}

template<class Type>
Foam::weightedMappedInternalValueFvPatchField<Type>::
weightedMappedInternalValueFvPatchField
(
    const weightedMappedInternalValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, fvMesh>& iF
)
:
    mappedInternalValueFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::weightedMappedInternalValueFvPatchField<Type>::updateCoeffs()
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

    const VolField<Type>& nbrField =
        this->mapper().sameRegion()
     && this->fieldName_ == this->internalField().name()
      ? refCast<const VolField<Type>>(this->internalField())
      : nbrMesh.template lookupObject<VolField<Type>>(this->fieldName_);

    // Construct mapped values
    Field<Type> sampleValues;

    if (this->interpolationScheme_ != interpolations::cell<Type>::typeName)
    {
        // Create an interpolation
        autoPtr<interpolation<Type>> interpolatorPtr
        (
            interpolation<Type>::New
            (
                this->interpolationScheme_,
                nbrField
            )
        );
        const interpolation<Type>& interpolator = interpolatorPtr();

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

    // Set the average, if necessary
    if (this->setAverage_)
    {
        const label patchi = this->patch().index();

        const surfaceScalarField::Boundary& phi =
            this->db().template lookupObject<surfaceScalarField>(phiName_)
           .boundaryField();

        const Type sampleAverageValue =
            gSum(phi[patchi]*this->patch().magSf()*sampleValues)
           /gSum(phi[patchi]*this->patch().magSf());

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

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::weightedMappedInternalValueFvPatchField<Type>::write(Ostream& os) const
{
    writeEntry(os, "phi", phiName_);

    mappedInternalValueFvPatchField<Type>::write(os);
}


// ************************************************************************* //

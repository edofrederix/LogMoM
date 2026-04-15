#ifndef weightedMappedInternalValueFvPatchField_H
#define weightedMappedInternalValueFvPatchField_H

#include "mappedInternalValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class weightedMappedInternalValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class weightedMappedInternalValueFvPatchField
:
    public mappedInternalValueFvPatchField<Type>
{
protected:

        //- Interpolation scheme to use
        const word phiName_;

public:

    //- Runtime type information
    TypeName("weightedMappedInternalValue");


    // Constructors

        //- Construct from patch, internal field and dictionary
        weightedMappedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given weightedMappedInternalValueFvPatchField
        //  onto a new patch
        weightedMappedInternalValueFvPatchField
        (
            const weightedMappedInternalValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        weightedMappedInternalValueFvPatchField
        (
            const weightedMappedInternalValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        weightedMappedInternalValueFvPatchField
        (
            const weightedMappedInternalValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new weightedMappedInternalValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "weightedMappedInternalValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

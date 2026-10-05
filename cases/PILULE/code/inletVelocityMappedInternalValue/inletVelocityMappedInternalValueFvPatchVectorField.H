#ifndef inletVelocityMappedInternalValueFvPatchVectorField_H
#define inletVelocityMappedInternalValueFvPatchVectorField_H

#include "mappedInternalValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class inletVelocityMappedInternalValueFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class inletVelocityMappedInternalValueFvPatchVectorField
:
    public mappedInternalValueFvPatchField<vector>
{
public:

    //- Runtime type information
    TypeName("inletVelocityMappedInternalValue");


    // Constructors

        //- Construct from patch, internal field and dictionary
        inletVelocityMappedInternalValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  inletVelocityMappedInternalValueFvPatchVectorField onto a new patch
        inletVelocityMappedInternalValueFvPatchVectorField
        (
            const inletVelocityMappedInternalValueFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        inletVelocityMappedInternalValueFvPatchVectorField
        (
            const inletVelocityMappedInternalValueFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        inletVelocityMappedInternalValueFvPatchVectorField
        (
            const inletVelocityMappedInternalValueFvPatchVectorField&,
            const DimensionedField<vector, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new inletVelocityMappedInternalValueFvPatchVectorField
                (
                    *this,
                    iF
                )
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

#endif

// ************************************************************************* //

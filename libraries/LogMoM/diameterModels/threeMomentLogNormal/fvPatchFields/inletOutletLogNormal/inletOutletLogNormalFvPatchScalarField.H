#ifndef inletOutletLogNormalFvPatchScalarField_H
#define inletOutletLogNormalFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class inletOutletLogNormalFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class inletOutletLogNormalFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private Data

        //- Value of desired sigma
        scalar sigma_;

        //- Value of desired Sauter mean diameter
        scalar dsm_;

        //- Name of alpha field
        word alphaName_;


protected:

    // Evaluation functions

            //- Compute specified moment
            tmp<scalarField> computeField() const;


public:

    //- Runtime type information
    TypeName("inletOutletLogNormal");


    // Constructors

        //- Construct from patch and internal field
        inletOutletLogNormalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        inletOutletLogNormalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  inletOutletLogNormalFvPatchScalarField onto a new patch
        inletOutletLogNormalFvPatchScalarField
        (
            const inletOutletLogNormalFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        inletOutletLogNormalFvPatchScalarField
        (
            const inletOutletLogNormalFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        inletOutletLogNormalFvPatchScalarField
        (
            const inletOutletLogNormalFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new inletOutletLogNormalFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


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

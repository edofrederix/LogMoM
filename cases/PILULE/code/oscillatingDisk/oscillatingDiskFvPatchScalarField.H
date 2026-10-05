#ifndef oscillatingDiskFvPatchScalarField_H
#define oscillatingDiskFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class oscillatingDiskFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class oscillatingDiskFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Spatial wavelength
        const scalar lambda_;

        //- Temporal wavelength
        const scalar kappa_;

        //- Average value
        const scalar average_;

        //- Wave amplitude
        const scalar amplitude_;

public:

    //- Runtime type information
    TypeName("oscillatingDisk");


    // Constructors

        //- Construct from patch, internal field and dictionary
        oscillatingDiskFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        oscillatingDiskFvPatchScalarField
        (
            const oscillatingDiskFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        oscillatingDiskFvPatchScalarField
        (
            const oscillatingDiskFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        oscillatingDiskFvPatchScalarField
        (
            const oscillatingDiskFvPatchScalarField&,
            const DimensionedField<scalar, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new oscillatingDiskFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Calculate and return the wavy field
            virtual tmp<scalarField> waves(const scalar t) const;

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

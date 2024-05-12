#ifndef BC_H
#define BC_H

#include "fixedValueFvPatchFields.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class TOPFLOWAlphaInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class TOPFLOWAlphaInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Average inlet value
        const scalar average_;

        //- Name of the phase (air or water)
        const word phaseName_;

        //- Case name
        const word caseName_;

        //- Raw data
        List<Tuple2<scalar,scalar>> data_;

public:

    //- Runtime type information
    TypeName("TOPFLOWAlphaInlet");


    // Constructors

        //- Construct from patch and internal field
        TOPFLOWAlphaInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        TOPFLOWAlphaInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        TOPFLOWAlphaInletFvPatchScalarField
        (
            const TOPFLOWAlphaInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        TOPFLOWAlphaInletFvPatchScalarField
        (
            const TOPFLOWAlphaInletFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        TOPFLOWAlphaInletFvPatchScalarField
        (
            const TOPFLOWAlphaInletFvPatchScalarField&,
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
                new TOPFLOWAlphaInletFvPatchScalarField
                (
                    *this,
                    iF
                )
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

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


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

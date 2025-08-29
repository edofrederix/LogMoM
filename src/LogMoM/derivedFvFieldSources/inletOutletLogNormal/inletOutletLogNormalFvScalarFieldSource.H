#ifndef inletOutletLogNormalFvScalarFieldSource_H
#define inletOutletLogNormalFvScalarFieldSource_H

#include "fvFieldSources.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class inletOutletLogNormalFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class inletOutletLogNormalFvScalarFieldSource
:
    public fvScalarFieldSource
{
private:

    // Private Data

        //- Value of desired sigma
        scalar sigma_;

        //- Value of desired Sauter-mean diameter
        scalar dsm_;

    // Private Member Functions

        //- Compute the q factor based on the moment dimension
        scalar q() const;


public:

    //- Runtime type information
    TypeName("inletOutletLogNormal");


    // Constructors

        //- Construct from internal field and dictionary
        inletOutletLogNormalFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        inletOutletLogNormalFvScalarFieldSource
        (
            const inletOutletLogNormalFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvFieldSource<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvFieldSource<scalar>>
            (
                new inletOutletLogNormalFvScalarFieldSource(*this, iF)
            );
        }


    //- Destructor
    virtual ~inletOutletLogNormalFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<DimensionedField<scalar, volMesh>> sourceValue
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the source value
        virtual tmp<Field<scalar>> sourceValue
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Return the source value
        virtual tmp<DimensionedField<scalar, volMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

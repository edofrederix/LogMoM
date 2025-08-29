#ifndef interfacialGrowthMomentFvScalarFieldSource_H
#define interfacialGrowthMomentFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class interfacialGrowthMomentFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class interfacialGrowthMomentFvScalarFieldSource
:
    public fvScalarFieldSource
{
protected:

    // Protected Member Functions

        //- Check and return the moment order of the field to which the field
        //  source is applied
        scalar gamma() const;

public:

    //- Runtime type information
    TypeName("interfacialGrowthMoment");


    // Constructors

        //- Construct from internal field and dictionary
        interfacialGrowthMomentFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        interfacialGrowthMomentFvScalarFieldSource
        (
            const interfacialGrowthMomentFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new interfacialGrowthMomentFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~interfacialGrowthMomentFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<DimensionedField<scalar, volMesh>> sourceValue
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the internal coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the source coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> sourceCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
            const DimensionedField<scalar, fvMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        interfacialGrowthMomentFvScalarFieldSource
        (
            const interfacialGrowthMomentFvScalarFieldSource&,
            const DimensionedField<scalar, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
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
        virtual tmp<DimensionedField<scalar, fvMesh>> sourceValue
        (
            const fvSource& model,
            const DimensionedField<scalar, fvMesh>& source
        ) const;

        //- Return the internal coefficient
        virtual tmp<DimensionedField<scalar, fvMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, fvMesh>& source
        ) const;

        //- Return the source coefficient
        virtual tmp<DimensionedField<scalar, fvMesh>> sourceCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, fvMesh>& source
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

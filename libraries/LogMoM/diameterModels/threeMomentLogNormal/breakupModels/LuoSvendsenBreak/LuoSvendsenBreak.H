#ifndef LuoSvendsenBreak_H
#define LuoSvendsenBreak_H

#include "breakupModel.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                       Class LuoSvendsenBreak Declaration
\*---------------------------------------------------------------------------*/

class LuoSvendsenBreak
:
    public breakupModel
{
private:

    // Private Data

        autoPtr<Function1s::Table<scalar>> gammaUpperReg5by11_;

        //- Empirical constant, defaults to 0.923
        const dimensionedScalar C4_;

        //- Empirical constant, defaults to 2.05
        const dimensionedScalar beta_;

        //- Ratio between minimum size of eddies in the inertial subrange
        //  and Kolmogorov length scale, defaults to 11.4
        const dimensionedScalar minEddyRatio_;

        //- Surface tension
        const dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("LuoSvendsen");

    // Constructor

        LuoSvendsenBreak
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LuoSvendsenBreak()
    {}

    // Member Functions

        // Return a field of binary breakup rates given two particle diameters
        virtual tmp<volScalarField> binaryRate
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

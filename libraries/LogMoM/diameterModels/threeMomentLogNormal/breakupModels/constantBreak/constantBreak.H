#ifndef constantBreak_H
#define constantBreak_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                        Class constantBreak Declaration
\*---------------------------------------------------------------------------*/

class constantBreak
:
    public breakupModel
{
    // Private Data

        //- Constant breakup rate
        dimensionedScalar B_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary
        constantBreak
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constantBreak();


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

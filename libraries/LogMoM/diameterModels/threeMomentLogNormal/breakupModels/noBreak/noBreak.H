#ifndef noBreak_H
#define noBreak_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                           Class noBreak Declaration
\*---------------------------------------------------------------------------*/

class noBreak
:
    public breakupModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary
        noBreak
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~noBreak();


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

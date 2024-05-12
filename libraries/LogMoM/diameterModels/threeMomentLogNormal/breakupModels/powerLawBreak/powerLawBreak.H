#ifndef powerLawBreak_H
#define powerLawBreak_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                         Class powerLawBreak Declaration
\*---------------------------------------------------------------------------*/

class powerLawBreak
:
    public breakupModel
{
    // Private Data

        //- Power
        dimensionedScalar power_;


public:

    //- Runtime type information
    TypeName("powerLaw");


    // Constructors

        //- Construct from a dictionary
        powerLawBreak
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~powerLawBreak();


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

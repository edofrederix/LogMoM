#ifndef powerLawUniform_H
#define powerLawUniform_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                         Class powerLawUniform Declaration
\*---------------------------------------------------------------------------*/

class powerLawUniform
:
    public breakupModel
{
    // Private Data

        //- Power
        dimensionedScalar power_;


public:

    //- Runtime type information
    TypeName("powerLawUniform");


    // Constructors

        //- Construct from a dictionary
        powerLawUniform
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~powerLawUniform();


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

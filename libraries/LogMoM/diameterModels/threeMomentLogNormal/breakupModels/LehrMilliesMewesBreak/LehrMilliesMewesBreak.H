#ifndef LehrMilliesMewesBreak_H
#define LehrMilliesMewesBreak_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                    Class LehrMilliesMewesBreak Declaration
\*---------------------------------------------------------------------------*/

class LehrMilliesMewesBreak
:
    public breakupModel
{
private:

    // Private Data

        //- Surface tension
        const dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("LehrMilliesMewes");

    // Constructor

        LehrMilliesMewesBreak
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LehrMilliesMewesBreak()
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

#ifndef LehrMilliesMewesCoaFreq_H
#define LehrMilliesMewesCoaFreq_H

#include "coalescenceFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceFrequencyModels
{

/*---------------------------------------------------------------------------*\
                           Class LehrMilliesMewesCoaFreq Declaration
\*---------------------------------------------------------------------------*/

class LehrMilliesMewesCoaFreq
:
    public coalescenceFrequencyModel
{
    // Private Data

        //- Critical velocity for coalescence, defaults to 0.08 m/s
        dimensionedScalar uCrit_;

        //- Maximum packing density of the bubbles, defaults to 0.6
        dimensionedScalar alphaMax_;

public:

    //- Runtime type information
    TypeName("LehrMilliesMewes");


    // Constructors

        //- Construct from a dictionary
        LehrMilliesMewesCoaFreq
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LehrMilliesMewesCoaFreq();


    // Member Functions

        // Return a field of coalescence frequencies given two particle sizes
        virtual tmp<volScalarField> frequency
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceFrequencyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

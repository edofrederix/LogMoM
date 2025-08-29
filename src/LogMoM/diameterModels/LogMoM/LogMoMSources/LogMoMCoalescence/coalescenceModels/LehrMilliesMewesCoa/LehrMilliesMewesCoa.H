#ifndef LehrMilliesMewes_H
#define LehrMilliesMewes_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                        Class LehrMilliesMewes Declaration
\*---------------------------------------------------------------------------*/

class LehrMilliesMewes
:
    public coalescenceModel
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
        LehrMilliesMewes
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LehrMilliesMewes();


    // Member Functions

        // Return a field of coalescence rates given two particle sizes
        virtual tmp<volScalarField> rate
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

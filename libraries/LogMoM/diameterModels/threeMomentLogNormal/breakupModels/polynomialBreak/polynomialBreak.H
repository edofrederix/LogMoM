#ifndef polynomialBreak_H
#define polynomialBreak_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                        Class polynomialBreak Declaration
\*---------------------------------------------------------------------------*/

class polynomialBreak
:
    public breakupModel
{
    // Private Data

        //- Polynomial coefficients
        List<scalar> B_;

        //- Polynomial powers
        List<scalar> r_;

public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from a dictionary
        polynomialBreak
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~polynomialBreak();


    // Member Functions

        // Return a field of breakup rates given a particle diameter
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

#ifndef polynomialCoaFreq_H
#define polynomialCoaFreq_H

#include "coalescenceFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceFrequencyModels
{

/*---------------------------------------------------------------------------*\
                     Class polynomialCoaFreq Declaration
\*---------------------------------------------------------------------------*/

class polynomialCoaFreq
:
    public coalescenceFrequencyModel
{
    // Private Data

        //- Polynomial coefficients
        List<scalar> K_;

        //- First polynomial power
        List<scalar> p_;

        //- Second polynomial power
        List<scalar> q_;


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from a dictionary
        polynomialCoaFreq
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~polynomialCoaFreq();


    // Member Functions

        // Return a field of coalescence frequencies given two particle sizes
        virtual tmp<volScalarField> frequency
        (
            const volScalarField& di,
            const volScalarField& dj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceFrequencyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

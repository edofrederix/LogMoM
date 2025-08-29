#ifndef polynomialCoa_H
#define polynomialCoa_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                        Class polynomialCoa Declaration
\*---------------------------------------------------------------------------*/

class polynomialCoa
:
    public coalescenceModel
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
        polynomialCoa
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~polynomialCoa();


    // Member Functions

        // Return a field of coalescence rates given two particle sizes
        virtual tmp<volScalarField> rate
        (
            const volScalarField& di,
            const volScalarField& dj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

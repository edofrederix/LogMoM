#ifndef polynomialCoaEff_H
#define polynomialCoaEff_H

#include "coalescenceEfficiencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceEfficiencyModels
{

/*---------------------------------------------------------------------------*\
                        Class polynomialCoaEff Declaration
\*---------------------------------------------------------------------------*/

class polynomialCoaEff
:
    public coalescenceEfficiencyModel
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
        polynomialCoaEff
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~polynomialCoaEff();


    // Member Functions

        // Return a field of coalescence efficiencies given two particle sizes
        virtual tmp<volScalarField> efficiency
        (
            const volScalarField& di,
            const volScalarField& dj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceEfficiencyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

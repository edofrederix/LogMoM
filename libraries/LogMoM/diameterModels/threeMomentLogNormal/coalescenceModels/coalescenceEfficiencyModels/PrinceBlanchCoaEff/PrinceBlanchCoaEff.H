#ifndef PrinceBlanchCoaEff_H
#define PrinceBlanchCoaEff_H

#include "coalescenceEfficiencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceEfficiencyModels
{

/*---------------------------------------------------------------------------*\
                      Class PrinceBlanchCoaEff Declaration
\*---------------------------------------------------------------------------*/

class PrinceBlanchCoaEff
:
    public coalescenceEfficiencyModel
{
    // Private Data

        //- Initial film thickness, defaults to 1e-4m
        dimensionedScalar h0_;

        //- Critical film thickness, defaults to 1e-8m
        dimensionedScalar hf_;

        //- Surface tension
        dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("PrinceBlanch");


    // Constructors

        //- Construct from a dictionary
        PrinceBlanchCoaEff
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~PrinceBlanchCoaEff();


    // Member Functions

        // Return a field of coalescence efficiencies given two particle sizes
        virtual tmp<volScalarField> efficiency
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceEfficiencyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

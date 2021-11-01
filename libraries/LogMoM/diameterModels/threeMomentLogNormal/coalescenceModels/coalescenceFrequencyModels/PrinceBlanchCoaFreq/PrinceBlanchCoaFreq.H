#ifndef PrinceBlanchCoaFreq_H
#define PrinceBlanchCoaFreq_H

#include "coalescenceFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceFrequencyModels
{

/*---------------------------------------------------------------------------*\
                           Class PrinceBlanchCoaFreq Declaration
\*---------------------------------------------------------------------------*/

class PrinceBlanchCoaFreq
:
    public coalescenceFrequencyModel
{
    // Private Data

        //- Optional coefficient C1, defaults to 0.356
        dimensionedScalar C1_;

        //- Surface tension
        dimensionedScalar sigma_;

        // Turbulence-driven coalescence
        Switch turbulent_;

        // Buoyancy-driven coalescence
        Switch buoyant_;

        // Laminar-shear-driven coalescence
        Switch laminar_;


public:

    //- Runtime type information
    TypeName("PrinceBlanch");


    // Constructors

        //- Construct from a dictionary
        PrinceBlanchCoaFreq
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~PrinceBlanchCoaFreq();


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

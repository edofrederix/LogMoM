#ifndef PrinceBlanch_H
#define PrinceBlanch_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                      Class PrinceBlanch Declaration
\*---------------------------------------------------------------------------*/

class PrinceBlanch
:
    public coalescenceModel
{
    // Private Data

        //- Initial film thickness, defaults to 1e-4m
        dimensionedScalar h0_;

        //- Critical film thickness, defaults to 1e-8m
        dimensionedScalar hf_;

        //- Optional coefficient C1, defaults to 0.356
        dimensionedScalar C1_;

        //- Switch for considering turbulent collisions
        Switch turbulence_;

        //- Switch for considering buoyancy-induced collisions
        Switch buoyancy_;

        //- Switch for considering buoyancy-induced collisions
        Switch laminarShear_;


public:

    //- Runtime type information
    TypeName("PrinceBlanch");


    // Constructors

        //- Construct from a dictionary
        PrinceBlanch
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~PrinceBlanch();


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

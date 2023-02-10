#ifndef constantCoaEff_H
#define constantCoaEff_H

#include "coalescenceEfficiencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceEfficiencyModels
{

/*---------------------------------------------------------------------------*\
                           Class constantCoaEff Declaration
\*---------------------------------------------------------------------------*/

class constantCoaEff
:
    public coalescenceEfficiencyModel
{
    // Private Data

        //- Constant coalescence efficiency
        dimensionedScalar K_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary
        constantCoaEff
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constantCoaEff();


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

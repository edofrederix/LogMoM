#ifndef noCoaEff_H
#define noCoaEff_H

#include "coalescenceEfficiencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceEfficiencyModels
{

/*---------------------------------------------------------------------------*\
                           Class noCoaEff Declaration
\*---------------------------------------------------------------------------*/

class noCoaEff
:
    public coalescenceEfficiencyModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary
        noCoaEff
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~noCoaEff();


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

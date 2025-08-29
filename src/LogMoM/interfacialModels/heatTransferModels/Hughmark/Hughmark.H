#ifndef Hughmark_H
#define Hughmark_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                        Class Hughmark Declaration
\*---------------------------------------------------------------------------*/

class Hughmark
:
    public heatTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    //- Runtime type information
    TypeName("Hughmark");


    // Constructors

        //- Construct from a dictionary and an interface
        Hughmark
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~Hughmark();


    // Member functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
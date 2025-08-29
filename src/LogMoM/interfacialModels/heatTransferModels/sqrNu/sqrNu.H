#ifndef sqrNu_H
#define sqrNu_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                        Class sqrNu Declaration
\*---------------------------------------------------------------------------*/

class sqrNu
:
    public heatTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;

        //- Nusselt value at zero point diameter
        dimensionedScalar Nu0_;

        //- Zero point diameter
        dimensionedScalar d0_;


public:

    //- Runtime type information
    TypeName("sqrNu");


    // Constructors

        //- Construct from a dictionary and an interface
        sqrNu
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~sqrNu();

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
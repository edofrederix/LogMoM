#ifndef NuHZDR_H
#define NuHZDR_H

#include "heatTransferModel.H"
#include "saturationTemperatureModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                        Class NuHZDR Declaration
\*---------------------------------------------------------------------------*/

class NuHZDR
:
    public heatTransferModel
{
    // Private Data

        //- Dispersed interface
        const dispersedPhaseInterface interface_;

        //- Saturation temperature model
        autoPtr<saturationTemperatureModel> saturationModelPtr_;

        //- Name of the mDot field
        word mDotName_;

        //- Name of the first phase of the phase change model, needed to
        //  adjust the sign of mDot if necessary
        word firstPhaseName_;


public:

    //- Runtime type information
    TypeName("NuHZDR");


    // Constructors

        //- Construct from a dictionary and an interface
        NuHZDR
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~NuHZDR();


    // Member functions

        // Return the continuous phase turbulence model
        virtual const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
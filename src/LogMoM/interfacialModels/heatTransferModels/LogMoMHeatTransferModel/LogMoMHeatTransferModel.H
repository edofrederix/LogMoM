#ifndef LogMoMHeatTransferModel_H
#define LogMoMHeatTransferModel_H

#include "heatTransferModel.H"
#include "LogMoMInterfacialModel.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                Class LogMoMHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMHeatTransferModel
:
    public heatTransferModel,
    public LogMoMInterfacialModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;

        //- The heat transfer model
        autoPtr<heatTransferModel> heatTransferModelPtr_;

        //- Residual alpha value used to pass argument
        scalar residualAlpha_;

        //- Assumed power of the Nusselt correlation, Nu ~ d^pow. Defaults to
        //  0.5.
        scalar pow_;

    // Private Member Functions

        //- Helper function to evaluate the scaled heat transfer function K of
        //  the underlying heat transfer model
        tmp<volScalarField> KdSqrPow() const;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMHeatTransferModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~LogMoMHeatTransferModel();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
#ifndef LogMoMTurbulentDispersionModel_H
#define LogMoMTurbulentDispersionModel_H

#include "dispersedTurbulentDispersionModel.H"
#include "LogMoMInterfacialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
                      Class LogMoMTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMTurbulentDispersionModel
:
    public dispersedTurbulentDispersionModel,
    public LogMoMInterfacialModel
{
    // Private Data

        //- The turbulent dispersion model
        autoPtr<turbulentDispersionModel> turbulentDispersionModel_;

    // Private Member Functions

        //- Turbulent diffusivity multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> DModel() const;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMTurbulentDispersionModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~LogMoMTurbulentDispersionModel();


    // Member Functions

        //- Return the turbulent diffusivity multiplying the gradient of the
        //  phase-fraction
        tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

#ifndef LogMoMPolyTurbulentDispersionModel_H
#define LogMoMPolyTurbulentDispersionModel_H

#include "polyTurbulentDispersionModel.H"
#include "LogMoMInterfacialModel.H"
#include "turbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
                      Class LogMoMPolyTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMPolyTurbulentDispersionModel
:
    public polyTurbulentDispersionModel,
    public LogMoMInterfacialModel
{
    //- Private Data

        //- The turbulent dispersion model
        autoPtr<turbulentDispersionModel> turbulentDispersionModel_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMPolyTurbulentDispersionModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~LogMoMPolyTurbulentDispersionModel();


    // Member Functions

        //- Turbulent diffusivity multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D() const;

        //- Turbulent diffusivity for moment gamma, multiplying the gradient of
        //  the phase-fraction
        virtual tmp<volScalarField> D
        (
            const label gamma
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

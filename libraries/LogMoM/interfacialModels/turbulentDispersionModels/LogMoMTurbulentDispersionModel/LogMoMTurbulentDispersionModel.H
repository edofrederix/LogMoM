#ifndef LogMoMTurbulentDispersionModel_H
#define LogMoMTurbulentDispersionModel_H

#include "turbulentDispersionModel.H"

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
    public turbulentDispersionModel
{
    // Private Data

        //- The turbulent dispersion model
        autoPtr<turbulentDispersionModel> turbulentDispersionModelPtr_;


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

        //- Turbulent diffusivity
        //  multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

#ifndef defaultPolyTurbulentDispersionModel_H
#define defaultPolyTurbulentDispersionModel_H

#include "polyTurbulentDispersionModel.H"
#include "turbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
            Class defaultPolyTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class defaultPolyTurbulentDispersionModel
:
    public polyTurbulentDispersionModel
{
    // Private Data

        //- The turbulent dispersion model
        autoPtr<turbulentDispersionModel> turbulentDispersionModelPtr_;


public:

    //- Runtime type information
    TypeName("default");


    // Constructors

        //- Construct from a dictionary and an interface
        defaultPolyTurbulentDispersionModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~defaultPolyTurbulentDispersionModel();


    // Member Functions

        //- Return the turbulent diffusivity multiplying the gradient of the
        //  phase-fraction
        tmp<volScalarField> D(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

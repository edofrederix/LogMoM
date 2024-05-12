#ifndef polyTurbulentDispersionModel_H
#define polyTurbulentDispersionModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class polyTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class polyTurbulentDispersionModel
{
public:

    //- Runtime type information
    TypeName("polyTurbulentDispersionModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            polyTurbulentDispersionModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );

    // Static Data Members

        //- Diffusivity dimensions
        static const dimensionSet dimD;

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = true;


    // Constructors

        //- Construct from a dictionary and an interface
        polyTurbulentDispersionModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~polyTurbulentDispersionModel();


    // Selectors

        static autoPtr<polyTurbulentDispersionModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );


    // Member Functions

        //- Return the turbulent diffusivity multiplying the gradient of the
        //  phase-fraction
        virtual tmp<volScalarField> D(const label gamma) const = 0;
};


/*---------------------------------------------------------------------------*\
                Class blendedPolyTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class blendedPolyTurbulentDispersionModel
:
    public BlendedInterfacialModel<polyTurbulentDispersionModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<polyTurbulentDispersionModel>::
            BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedPolyTurbulentDispersionModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the turbulent diffusivity multiplying the gradient of the
        //  phase-fraction
        tmp<volScalarField> D(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

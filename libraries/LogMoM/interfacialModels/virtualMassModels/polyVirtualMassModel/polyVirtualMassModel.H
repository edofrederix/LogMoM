#ifndef polyVirtualMassModel_H
#define polyVirtualMassModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class polyVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class polyVirtualMassModel
:
    public regIOobject
{
public:

    //- Runtime type information
    TypeName("polyVirtualMassModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            polyVirtualMassModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface,
                const bool registerObject
            ),
            (dict, interface, registerObject)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimK;

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = true;


    // Constructors

        //- Construct from a dictionary and an interface
        polyVirtualMassModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~polyVirtualMassModel();


    // Selectors

        static autoPtr<polyVirtualMassModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true,
            const bool registerObject=false
        );


    // Member Functions

        //- Return the virtual mass coefficient K
        //  used in the momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... K*(DU1_Dt - DU2_Dt)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(DU1_Dt - DU2_Dt)
        virtual tmp<volScalarField> K(const label gamma) const = 0;

        //- Return the virtual mass coefficient Kf
        //  used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf(const label gamma) const = 0;

        // Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                      Class blendedPolyVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class blendedPolyVirtualMassModel
:
    public BlendedInterfacialModel<polyVirtualMassModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<polyVirtualMassModel>
          ::BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedPolyVirtualMassModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the lift coefficient K
        tmp<volScalarField> K(const label gamma) const;

        //- Return the lift coefficient Kf
        tmp<surfaceScalarField> Kf(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

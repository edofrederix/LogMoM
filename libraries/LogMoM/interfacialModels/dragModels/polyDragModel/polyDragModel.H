#ifndef polyDragModel_H
#define polyDragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

namespace Foam
{

class swarmCorrection;

/*---------------------------------------------------------------------------*\
                          Class polyDragModel Declaration
\*---------------------------------------------------------------------------*/

class polyDragModel
:
    public regIOobject
{
public:

    //- Runtime type information
    TypeName("polyDragModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            polyDragModel,
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

        // Construct from a dictionary and an interface
        polyDragModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~polyDragModel();


    // Selectors

        static autoPtr<polyDragModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true,
            const bool registerObject=false
        );


    // Member Functions

        //- Return the drag coefficient K
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(U2-U1)
        virtual tmp<volScalarField> K
        (
            const label gamma
        ) const = 0;

        //- Return the drag coefficient Kf
        //  used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf
        (
            const label gamma
        ) const = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                      Class blendedPolyDragModel Declaration
\*---------------------------------------------------------------------------*/

class blendedPolyDragModel
:
    public BlendedInterfacialModel<polyDragModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using BlendedInterfacialModel<polyDragModel>::BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedPolyDragModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the drag coefficient K
        tmp<volScalarField> K
        (
            const label gamma
        ) const;

        //- Return the drag coefficient Kf
        tmp<surfaceScalarField> Kf
        (
            const label gamma
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

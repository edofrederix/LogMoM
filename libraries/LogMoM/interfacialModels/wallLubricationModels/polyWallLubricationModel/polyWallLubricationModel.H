#ifndef polyWallLubricationModel_H
#define polyWallLubricationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallDependentModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class polyWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class polyWallLubricationModel
:
    public wallDependentModel
{
public:

    //- Runtime type information
    TypeName("polyWallLubricationModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            polyWallLubricationModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimF;

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = true;


    // Constructors

        //- Construct from a dictionary and an interface
        polyWallLubricationModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~polyWallLubricationModel();


    // Selectors

        static autoPtr<polyWallLubricationModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );


    // Member Functions

        //- Return wall lubrication force
        virtual tmp<volVectorField> F(const label gamma) const = 0;

        //- Return face wall lubrication force
        virtual tmp<surfaceScalarField> Ff(const label gamma) const = 0;
};


/*---------------------------------------------------------------------------*\
                   Class blendedPolyWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class blendedPolyWallLubricationModel
:
    public BlendedInterfacialModel<polyWallLubricationModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<polyWallLubricationModel>::
            BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedPolyWallLubricationModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return wall lubrication force
        tmp<volVectorField> F(const label gamma) const;

        //- Return face wall lubrication force
        tmp<surfaceScalarField> Ff(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

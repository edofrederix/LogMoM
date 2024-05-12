#ifndef polyLiftModel_H
#define polyLiftModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class polyLiftModel Declaration
\*---------------------------------------------------------------------------*/

class polyLiftModel
{
public:

    //- Runtime type information
    TypeName("polyLiftModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            polyLiftModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Static Data Members

        //- Force dimensions
        static const dimensionSet dimF;

        //- Does this model require correcting on fixed flux boundaries?
        static const bool correctFixedFluxBCs = true;


    // Constructors

        //- Construct from a dictionary and an interface
        polyLiftModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~polyLiftModel();


    // Selectors

        static autoPtr<polyLiftModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );


    // Member Functions

        //- Return lift force
        virtual tmp<volVectorField> F(const label gamma) const = 0;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff(const label gamma) const = 0;
};


/*---------------------------------------------------------------------------*\
                      Class blendedPolyLiftModel Declaration
\*---------------------------------------------------------------------------*/

class blendedPolyLiftModel
:
    public BlendedInterfacialModel<polyLiftModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using BlendedInterfacialModel<polyLiftModel>::BlendedInterfacialModel;


    // Selectors

        static autoPtr<blendedPolyLiftModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return lift force
        tmp<volVectorField> F(const label gamma) const;

        //- Return face lift force
        tmp<surfaceScalarField> Ff(const label gamma) const;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

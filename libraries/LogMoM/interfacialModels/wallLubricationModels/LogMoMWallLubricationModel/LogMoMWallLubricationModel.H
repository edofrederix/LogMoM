#ifndef LogMoMWallLubricationModel_H
#define LogMoMWallLubricationModel_H

#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class LogMoMWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMWallLubricationModel
:
    public wallLubricationModel
{
    // Private Data

        //- The wall lubrication model
        autoPtr<wallLubricationModel> wallLubricationModelPtr_;


public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMWallLubricationModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~LogMoMWallLubricationModel();


    // Member Functions

        //- Return wall lubrication force
        virtual tmp<volVectorField> F() const;

        //- Return face wall lubrication force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

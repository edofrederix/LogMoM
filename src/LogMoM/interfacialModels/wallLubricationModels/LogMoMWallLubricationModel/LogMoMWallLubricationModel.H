#ifndef LogMoMWallLubricationModel_H
#define LogMoMWallLubricationModel_H

#include "dispersedWallLubricationModel.H"
#include "LogMoMInterfacialModel.H"
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
    public dispersedWallLubricationModel,
    public LogMoMInterfacialModel
{
    // Private Data

        //- The wall lubrication model
        autoPtr<wallLubricationModel> wallLubricationModel_;

    // Private Member Functions

        //- Return phase-intensive wall lubrication force of the actual
        //  model
        tmp<volVectorField> FiModel() const;

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

        //- Return phase-intensive wall lubrication force. Note that this is
        //  signed relative to dispersed-continuous, not phase1-phase2.
        virtual tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

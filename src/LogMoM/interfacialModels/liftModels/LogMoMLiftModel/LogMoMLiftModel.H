#ifndef LogMoMLiftModel_H
#define LogMoMLiftModel_H

#include "dispersedLiftModel.H"
#include "LogMoMInterfacialModel.H"
#include "wallDampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                    Class LogMoMLiftModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMLiftModel
:
    public dispersedLiftModel,
    public LogMoMInterfacialModel
{
private:

    // Private Data

        //- The lift model
        autoPtr<liftModel> liftModelPtr_;

        //- Optional wall-damping model
        autoPtr<wallDampingModel> wallDampingModelPtr_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMLiftModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~LogMoMLiftModel();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;

        //- Return phase-intensive lift force. Note that this is
        //  signed relative to dispersed-continuous, not phase1-phase2.
        virtual tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

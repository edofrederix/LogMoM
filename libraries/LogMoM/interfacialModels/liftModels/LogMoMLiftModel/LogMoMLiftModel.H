#ifndef LogMoMLiftModel_H
#define LogMoMLiftModel_H

#include "liftModel.H"
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
    public liftModel
{
private:

    //- Private Data

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

        //- Return lift force
        virtual tmp<volVectorField> F() const;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

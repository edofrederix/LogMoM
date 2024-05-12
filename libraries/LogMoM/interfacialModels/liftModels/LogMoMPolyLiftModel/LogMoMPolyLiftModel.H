#ifndef LogMoMPolyLiftModel_H
#define LogMoMPolyLiftModel_H

#include "polyLiftModel.H"
#include "LogMoMInterfacialModel.H"
#include "liftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                    Class LogMoMPolyLiftModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMPolyLiftModel
:
    public polyLiftModel,
    public LogMoMInterfacialModel
{
private:

    //- Private Data

        //- The lift model
        autoPtr<liftModel> liftModelPtr_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMPolyLiftModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~LogMoMPolyLiftModel();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;

        //- Lift coefficient multiplied by Ur
        virtual tmp<volVectorField> ClUr() const;

        //- Return phase-intensive lift force
        virtual tmp<volVectorField> Fi() const;

        //- Return phase-intensive lift force for moment gamma
        virtual tmp<volVectorField> Fi
        (
            const label gamma
        ) const;

        //- Return lift force
        virtual tmp<volVectorField> F(const label gamma) const;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

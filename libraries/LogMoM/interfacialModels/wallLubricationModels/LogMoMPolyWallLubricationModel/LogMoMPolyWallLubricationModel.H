#ifndef LogMoMPolyWallLubricationModel_H
#define LogMoMPolyWallLubricationModel_H

#include "polyWallLubricationModel.H"
#include "LogMoMInterfacialModel.H"
#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                     Class LogMoMPolyWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMPolyWallLubricationModel
:
    public polyWallLubricationModel,
    public LogMoMInterfacialModel
{
    // Private Data

        //- The wall lubrication model
        autoPtr<wallLubricationModel> wallLubricationModel_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMPolyWallLubricationModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~LogMoMPolyWallLubricationModel();


    // Member Functions

        //- Wall lubrication force from the underlying model
        tmp<volVectorField> Fim() const;

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi() const;

        //- Return moment-intensive wall lubrication force
        tmp<volVectorField> Fi
        (
            const label gamma
        ) const;

        //- Return wall lubrication force
        virtual tmp<volVectorField> F(const label gamma) const;

        //- Return face wall lubrication force
        virtual tmp<surfaceScalarField> Ff(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

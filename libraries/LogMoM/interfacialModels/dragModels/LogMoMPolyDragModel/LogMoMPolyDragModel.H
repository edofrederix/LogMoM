#ifndef LogMoMPolyDragModel_H
#define LogMoMPolyDragModel_H

#include "polyDragModel.H"
#include "LogMoMInterfacialModel.H"
#include "dragModel.H"
#include "swarmCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                       Class LogMoMPolyDragModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMPolyDragModel
:
    public polyDragModel,
    public LogMoMInterfacialModel
{
    // Private Data

        //- The drag model
        autoPtr<dragModel> dragModelPtr_;

        //- Swarm correction
        autoPtr<swarmCorrection> swarmCorrection_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMPolyDragModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~LogMoMPolyDragModel();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;

        //- Drag coefficient scaled by d^2
        virtual tmp<volScalarField> CdReBydSqr() const;

        //- Return the phase-intensive drag coefficient Ki used in the momentum
        //  equations
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(U2-U1)
        virtual tmp<volScalarField> Ki() const;

        //- Return the moment-intensive drag coefficient Ki used in the momentum
        //  equations for moment gamma
        virtual tmp<volScalarField> Ki
        (
            const label gamma
        ) const;

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

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

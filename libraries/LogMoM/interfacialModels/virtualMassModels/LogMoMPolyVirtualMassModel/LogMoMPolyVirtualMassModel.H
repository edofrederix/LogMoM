#ifndef LogMoMPolyVirtualMassModel_H
#define LogMoMPolyVirtualMassModel_H

#include "polyVirtualMassModel.H"
#include "LogMoMInterfacialModel.H"
#include "virtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace virtualMassModels
{

/*---------------------------------------------------------------------------*\
                        Class LogMoMPolyVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMPolyVirtualMassModel
:
    public polyVirtualMassModel,
    public LogMoMInterfacialModel
{
private:

    //- Private Data

        //- The virtual mass model
        autoPtr<virtualMassModel> virtualMassModelPtr_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMPolyVirtualMassModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~LogMoMPolyVirtualMassModel();


    // Member Functions

        //- Virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const;

        //- Return the phase-intensive virtual mass coefficient Ki used in the
        //  momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(DU1_Dt - DU2_Dt)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(DU1_Dt - DU2_Dt)
        virtual tmp<volScalarField> Ki() const;

        //- Return the moment-intensive virtual mass coefficient Ki used in the
        //  momentum equation for moment gamma
        virtual tmp<volScalarField> Ki
        (
            const label gamma
        ) const;

        //- Return the lift coefficient K
        tmp<volScalarField> K(const label gamma) const;

        //- Return the lift coefficient Kf
        tmp<surfaceScalarField> Kf(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace virtualMassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

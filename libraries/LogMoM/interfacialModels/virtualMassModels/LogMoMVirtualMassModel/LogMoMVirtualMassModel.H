#ifndef LogMoMVirtualMassModel_H
#define LogMoMVirtualMassModel_H

#include "virtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace virtualMassModels
{

/*---------------------------------------------------------------------------*\
                            Class LogMoMVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMVirtualMassModel
:
    public virtualMassModel
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
        LogMoMVirtualMassModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~LogMoMVirtualMassModel();


    // Member Functions

        //- Return the virtual mass coefficient K
        //  used in the momentum equation
        //    ddt(alpha1*rho1*U1) + ... = ... K*(DU1_Dt - DU2_Dt)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(DU1_Dt - DU2_Dt)
        virtual tmp<volScalarField> K() const;

        //- Return the virtual mass coefficient Kf
        //  used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace virtualMassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

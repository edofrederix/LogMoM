#ifndef LogMoMVirtualMassModel_H
#define LogMoMVirtualMassModel_H

#include "dispersedVirtualMassModel.H"
#include "LogMoMInterfacialModel.H"
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
    public dispersedVirtualMassModel,
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
        LogMoMVirtualMassModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~LogMoMVirtualMassModel();


    // Member Functions

        //- Return the virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const;

        //- Return the phase-intensive virtual mass coefficient Ki
        virtual tmp<volScalarField> Ki() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace virtualMassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

#ifndef defaultPolyVirtualMassModel_H
#define defaultPolyVirtualMassModel_H

#include "polyVirtualMassModel.H"
#include "virtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace virtualMassModels
{

/*---------------------------------------------------------------------------*\
                  Class defaultPolyVirtualMassModel Declaration
\*---------------------------------------------------------------------------*/

class defaultPolyVirtualMassModel
:
    public polyVirtualMassModel
{
private:

    //- Private Data

        //- The virtual mass model
        autoPtr<virtualMassModel> virtualMassModelPtr_;


public:

    //- Runtime type information
    TypeName("default");


    // Constructors

        //- Construct from a dictionary and an interface
        defaultPolyVirtualMassModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~defaultPolyVirtualMassModel();


    // Member Functions

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

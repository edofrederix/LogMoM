#ifndef LogMoMDragModel_H
#define LogMoMDragModel_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                       Class LogMoMDragModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMDragModel
:
    public dragModel
{
    // Private Data

        //- The drag model
        autoPtr<dragModel> dragModelPtr_;


public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMDragModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~LogMoMDragModel();


    // Member Functions

        //- Return the drag coefficient K
        tmp<volScalarField> K() const;

        //- Return the drag coefficient Kf
        tmp<surfaceScalarField> Kf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

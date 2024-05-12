#ifndef defaultPolyDragModel_H
#define defaultPolyDragModel_H

#include "polyDragModel.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                       Class defaultPolyDragModel Declaration
\*---------------------------------------------------------------------------*/

class defaultPolyDragModel
:
    public polyDragModel
{
    // Private Data

        //- Actual drag model
        autoPtr<dragModel> dragModelPtr_;


public:

    //- Runtime type information
    TypeName("default");


    // Constructors

        //- Construct from a dictionary and an interface
        defaultPolyDragModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~defaultPolyDragModel();


    // Member Functions

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

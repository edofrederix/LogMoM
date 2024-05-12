#ifndef defaultPolyWallLubricationModel_H
#define defaultPolyWallLubricationModel_H

#include "polyWallLubricationModel.H"
#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class defaultPolyWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class defaultPolyWallLubricationModel
:
    public polyWallLubricationModel
{
    // Private Data

        //- The wall lubrication model
        autoPtr<wallLubricationModel> wallLubricationModelPtr_;


public:

    //- Runtime type information
    TypeName("default");


    // Constructors

        //- Construct from a dictionary and an interface
        defaultPolyWallLubricationModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~defaultPolyWallLubricationModel();


        //- Return wall lubrication force
        tmp<volVectorField> F(const label gamma) const;

        //- Return face wall lubrication force
        tmp<surfaceScalarField> Ff(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

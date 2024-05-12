#ifndef defaultPolyLiftModel_H
#define defaultPolyLiftModel_H

#include "polyLiftModel.H"
#include "liftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                        Class defaultPolyLiftModel Declaration
\*---------------------------------------------------------------------------*/

class defaultPolyLiftModel
:
    public polyLiftModel
{
private:

    //- Private Data

        //- The lift model
        autoPtr<liftModel> liftModelPtr_;


public:

    //- Runtime type information
    TypeName("default");


    // Constructors

        //- Construct from a dictionary and an interface
        defaultPolyLiftModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~defaultPolyLiftModel();


    // Member Functions

        //- Return lift force
        tmp<volVectorField> F(const label gamma) const;

        //- Return face lift force
        tmp<surfaceScalarField> Ff(const label gamma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

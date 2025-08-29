#ifndef SugrueLift_H
#define SugrueLift_H

#include "dispersedLiftModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                        Class SugrueLift Declaration
\*---------------------------------------------------------------------------*/

class SugrueLift
:
    public dispersedLiftModel
{
private:

    //- Private Data

        //- The aspect ratio model
        autoPtr<aspectRatioModel> aspectRatio_;

public:

    //- Runtime type information
    TypeName("Sugrue");


    // Constructors

        //- Construct from a dictionary and an interface
        SugrueLift
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~SugrueLift();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;

        //- Return a reference to the turbulence model for the continuous phase
        const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

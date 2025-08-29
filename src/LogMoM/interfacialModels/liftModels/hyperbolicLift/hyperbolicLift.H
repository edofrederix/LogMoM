#ifndef hyperbolicLift_H
#define hyperbolicLift_H

#include "dispersedLiftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aspectRatioModel;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                        Class hyperbolicLift Declaration
\*---------------------------------------------------------------------------*/

class hyperbolicLift
:
    public dispersedLiftModel
{
private:

    //- Private Data

        //- Lower and upper lift coefficients
        scalar C1_;
        scalar C2_;

        //- Step diameter
        dimensionedScalar d0_;

        //- Step slope
        scalar C_;


public:

    //- Runtime type information
    TypeName("hyperbolic");


    // Constructors

        //- Construct from a dictionary and an interface
        hyperbolicLift
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~hyperbolicLift();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

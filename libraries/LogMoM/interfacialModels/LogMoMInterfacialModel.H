#ifndef LogMoMInterfacialModel_H
#define LogMoMInterfacialModel_H

#include "dispersedDragModel.H"
#include "GaussQuadrature.H"
#include "threeMomentLogNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class LogMoMInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

class LogMoMInterfacialModel
{
public:

    // Public data types

        //- Integration mode enumeration
        enum class integrationMode
        {
            quadrature,
            Riemann,
            diameter
        };

        //- Integration mode names
        static const NamedEnum<integrationMode, 3> integrationModeNames_;

private:

    // Private Data

        //- Integration mode
        integrationMode integrationMode_;

        //- For the diameter integration mode we need the p and q coefficients
        //  of the average diameter d_{p,q}, which are provided as polynomials
        //  of gamma with the first entry multiplying gamma^0, the second entry
        //  gamma^1, etc.
        scalarList p_;
        scalarList q_;

        //- By default we have p = {0,1} and q = {2,0}, i.e., d_{gamma,2}, so
        //  that the Sauter mean diameter is obtained for gamma = 3.
        static const scalarList pd_;
        static const scalarList qd_;

        //- Number of Riemann bins
        label Nr_;

        //- Blend velocity
        Switch blendVelocity_;

protected:

    // Protected Data

        //- The dispersed phase interface
        dispersedPhaseInterface interface_;

        //- Gauss-Hermite quadrature rules
        autoPtr<GaussQuadrature> quadrature_;

        //- The LogMoM diameter model
        diameterModels::threeMomentLogNormal* diameterModelPtr_;

public:

    //- Runtime type information
    TypeName("LogMoM");


    // Constructors

        //- Construct from a dictionary and an interface
        LogMoMInterfacialModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~LogMoMInterfacialModel();


    // Member Functions

        //- Evaluate the function (1/M_p) int [d^p*f()*n(d)] dd on [0,inf],
        //  where f() is provided as a function pointer and p is a polynomial
        //  function of gamma. If the integration mode is 'quadrature',
        //  Gauss-Hermite quadrature is used. If the integration mode is
        //  'diameter', a specified mean diameter d_{p,q} is used as an
        //  integration point to evaluate f().
        template<class Type, class Model>
        tmp<GeometricField<Type, fvPatchField, volMesh>> evaluate
        (
            const label gamma,
            tmp<GeometricField<Type, fvPatchField, volMesh>>
                (Model::*f)() const,
            const Model& model
        ) const;

        // Return a copy of the dispersed phase velocity field of moment gamma
        tmp<volVectorField> U(const label gamma) const;

        // Linear blend of the three moment velocities defined at their
        // respective diameters. When the diameter is outside the range [d0,
        // d3], the velocity is limited to U0 or U3.
        tmp<volVectorField> UBlend
        (
            const volScalarField& d,
            const volScalarField& d0,
            const volScalarField& d2,
            const volScalarField& d3,
            const volVectorField& U0,
            const volVectorField& U2,
            const volVectorField& U3
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LogMoMInterfacialModelTemplates.C"

#endif

// ************************************************************************* //

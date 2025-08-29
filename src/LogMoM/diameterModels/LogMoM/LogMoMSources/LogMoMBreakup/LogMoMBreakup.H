#ifndef LogMoMBreakup_H
#define LogMoMBreakup_H

#include "LogMoMSource.H"
#include "breakupModel.H"
#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace LogMoMSources
{

/*---------------------------------------------------------------------------*\
                      Class LogMoMBreakup Declaration
\*---------------------------------------------------------------------------*/

class LogMoMBreakup
:
    public LogMoMSource
{
private:

    // Private Data

        // Gauss-Hermite quadrature
        autoPtr<GaussQuadrature> GHQ_;

        // Gauss-Legendre quadrature
        autoPtr<GaussQuadrature> GLQ_;

        //- Breakup model
        autoPtr<breakupModel> breakupModel_;


public:

    //- Runtime type information
    TypeName("breakup");

    // Constructors

        LogMoMBreakup
        (
            const LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LogMoMBreakup()
    {}


    // Member Functions

        // Moment source term
        virtual tmp<fvScalarMatrix> R(volScalarField&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LogMoMSources
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

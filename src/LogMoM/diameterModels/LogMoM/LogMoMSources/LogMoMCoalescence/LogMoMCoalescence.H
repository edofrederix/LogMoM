#ifndef LogMoMCoalescence_H
#define LogMoMCoalescence_H

#include "LogMoMSource.H"
#include "coalescenceModel.H"
#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace LogMoMSources
{

/*---------------------------------------------------------------------------*\
                      Class LogMoMCoalescence Declaration
\*---------------------------------------------------------------------------*/

class LogMoMCoalescence
:
    public LogMoMSource
{
private:

    // Private Data

        // Gauss-Hermite quadrature
        autoPtr<GaussQuadrature> GHQ_;

        //- Coalescence model
        autoPtr<coalescenceModel> coalescenceModel_;

public:

    //- Runtime type information
    TypeName("coalescence");

    // Constructors

        LogMoMCoalescence
        (
            const LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LogMoMCoalescence()
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

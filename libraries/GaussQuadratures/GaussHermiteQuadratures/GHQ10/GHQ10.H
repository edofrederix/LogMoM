#ifndef GHQ10_H
#define GHQ10_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                            Class GHQ10 Declaration
\*---------------------------------------------------------------------------*/

class GHQ10
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GHQ10");


    // Constructors
    GHQ10();


    //- Destructor
    virtual ~GHQ10();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadratures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

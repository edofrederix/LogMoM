#ifndef GHQ20_H
#define GHQ20_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                           Class GHQ20 Declaration
\*---------------------------------------------------------------------------*/

class GHQ20
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GHQ20");


    // Constructors
    GHQ20();


    //- Destructor
    virtual ~GHQ20();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadratures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

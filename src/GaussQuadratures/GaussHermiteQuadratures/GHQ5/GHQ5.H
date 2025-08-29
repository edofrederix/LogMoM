#ifndef GHQ5_H
#define GHQ5_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                            Class GHQ5 Declaration
\*---------------------------------------------------------------------------*/

class GHQ5
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GHQ5");


    // Constructors
    GHQ5();


    //- Destructor
    virtual ~GHQ5();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadrature
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

#ifndef GHQ3_H
#define GHQ3_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                             Class GHQ5 Declaration
\*---------------------------------------------------------------------------*/

class GHQ3
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GHQ3");


    // Constructors
    GHQ3();


    //- Destructor
    virtual ~GHQ3();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadrature
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

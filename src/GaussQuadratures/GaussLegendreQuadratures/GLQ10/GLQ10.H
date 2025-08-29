#ifndef GLQ10_H
#define GLQ10_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                            Class GLQ10 Declaration
\*---------------------------------------------------------------------------*/

class GLQ10
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GLQ10");


    // Constructors
    GLQ10();


    //- Destructor
    virtual ~GLQ10();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadratures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

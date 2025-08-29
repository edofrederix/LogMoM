#ifndef GLQ20_H
#define GLQ20_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                             Class GLQ20 Declaration
\*---------------------------------------------------------------------------*/

class GLQ20
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GLQ20");


    // Constructors
    GLQ20();


    //- Destructor
    virtual ~GLQ20();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadratures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

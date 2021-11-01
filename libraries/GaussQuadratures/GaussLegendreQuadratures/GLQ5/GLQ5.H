#ifndef GLQ5_H
#define GLQ5_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                            Class GLQ5 Declaration
\*---------------------------------------------------------------------------*/

class GLQ5
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GLQ5");


    // Constructors
    GLQ5();


    //- Destructor
    virtual ~GLQ5();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadratures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

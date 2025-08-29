#ifndef GLQ3_H
#define GLQ3_H

#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace GaussQuadratures
{

/*---------------------------------------------------------------------------*\
                             Class GLQ3 Declaration
\*---------------------------------------------------------------------------*/

class GLQ3
:
    public GaussQuadrature
{

public:

    //- Runtime type information
    TypeName("GLQ3");


    // Constructors
    GLQ3();


    //- Destructor
    virtual ~GLQ3();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GaussQuadratures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

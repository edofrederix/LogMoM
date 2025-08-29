#include "GLQ3.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace GaussQuadratures
{
    defineTypeNameAndDebug(GLQ3, 0);
    addToRunTimeSelectionTable(GaussQuadrature, GLQ3, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GLQ3::GLQ3()
:
    GaussQuadrature()
{
    x_.setSize(3);
    w_.setSize(3);

    x_[0] = -0.7745966692;
    x_[1] = 0.0;
    x_[2] = 0.7745966692;

    w_[0] = 0.5555555556;
    w_[1] = 0.8888888889;
    w_[2] = 0.5555555556;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GLQ3::~GLQ3()
{}

// ************************************************************************* //

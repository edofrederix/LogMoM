#include "GLQ5.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace GaussQuadratures
{
    defineTypeNameAndDebug(GLQ5, 0);
    addToRunTimeSelectionTable(GaussQuadrature, GLQ5, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GLQ5::GLQ5()
:
    GaussQuadrature()
{
    x_.setSize(5);
    w_.setSize(5);

    x_[0] = -0.9061798459386640;
    x_[1] = -0.5384693101056831;
    x_[2] = 0.0;
    x_[3] = 0.5384693101056831;
    x_[4] = 0.9061798459386640;

    w_[0] = 0.2369268850561891;
    w_[1] = 0.4786286704993665;
    w_[2] = 0.5688888888888889;
    w_[3] = 0.4786286704993665;
    w_[4] = 0.2369268850561891;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GLQ5::~GLQ5()
{}

// ************************************************************************* //

#include "GHQ5.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace GaussQuadratures
{
    defineTypeNameAndDebug(GHQ5, 0);
    addToRunTimeSelectionTable(GaussQuadrature, GHQ5, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GHQ5::GHQ5()
:
    GaussQuadrature()
{
    x_.setSize(5);
    w_.setSize(5);

    x_[0] = -0.202018287045609e1;
    x_[1] = -0.958572464613819;
    x_[2] =  0.0;
    x_[3] =  0.958572464613819;
    x_[4] =  0.202018287045609e1;

    w_[0] =  0.199532420590459e-1;
    w_[1] =  0.393619323152241;
    w_[2] =  0.945308720482942;
    w_[3] =  0.393619323152241;
    w_[4] =  0.199532420590459e-1;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GHQ5::~GHQ5()
{}

// ************************************************************************* //

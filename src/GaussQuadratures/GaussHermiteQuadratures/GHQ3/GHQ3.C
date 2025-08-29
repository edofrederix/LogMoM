#include "GHQ3.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace GaussQuadratures
{
    defineTypeNameAndDebug(GHQ3, 0);
    addToRunTimeSelectionTable(GaussQuadrature, GHQ3, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GHQ3::GHQ3()
:
    GaussQuadrature()
{
    x_.setSize(3);
    w_.setSize(3);

    x_[0] = -1.224744871391589049099;
    x_[1] =  0.0;
    x_[2] =  1.224744871391589049099;

    w_[0] =  0.295408975150919337883;
    w_[1] =  1.181635900603677351532;
    w_[2] =  0.295408975150919337883;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GaussQuadratures::GHQ3::~GHQ3()
{}

// ************************************************************************* //

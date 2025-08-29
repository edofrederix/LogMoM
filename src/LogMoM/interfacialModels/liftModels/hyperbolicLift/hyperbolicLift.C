#include "hyperbolicLift.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(hyperbolicLift, 0);
    addToRunTimeSelectionTable(liftModel, hyperbolicLift, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::hyperbolicLift::hyperbolicLift
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedLiftModel(dict, interface),
    C1_(readScalar(dict.lookup("C1"))),
    C2_(readScalar(dict.lookup("C2"))),
    d0_(dimLength, readScalar(dict.lookup("d"))),
    C_(readScalar(dict.lookup("C")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::hyperbolicLift::~hyperbolicLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::hyperbolicLift::Cl() const
{
    const volScalarField& d = interface_.dispersed().d();

    return (C2_-C1_)/2.0 * tanh(C_*log(d/d0_)) + (C1_+C2_)/2.0;
}


// ************************************************************************* //

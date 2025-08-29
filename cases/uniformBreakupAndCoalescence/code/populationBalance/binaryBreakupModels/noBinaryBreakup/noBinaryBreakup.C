#include "noBinaryBreakup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(noBinaryBreakup, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        noBinaryBreakup,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::binaryBreakupModels::noBinaryBreakup::
noBinaryBreakup
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::binaryBreakupModels::noBinaryBreakup::
addToBinaryBreakupRate
(
    volScalarField::Internal& binaryBreakupRate,
    const label i,
    const label j
)
{}


// ************************************************************************* //

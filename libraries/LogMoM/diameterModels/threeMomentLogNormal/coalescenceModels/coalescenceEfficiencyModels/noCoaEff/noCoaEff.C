#include "noCoaEff.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceEfficiencyModels
{
    defineTypeNameAndDebug(noCoaEff, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyModel,
        noCoaEff,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::noCoaEff::noCoaEff
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    coalescenceEfficiencyModel(pair, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::noCoaEff::~noCoaEff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyModels::noCoaEff::efficiency
(
    const volScalarField& di,
    const volScalarField& dj
) const
{
    tmp<volScalarField> tE
    (
        volScalarField::New
        (
            "tE",
            pair_.mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );

    return tE;
}

// ************************************************************************* //

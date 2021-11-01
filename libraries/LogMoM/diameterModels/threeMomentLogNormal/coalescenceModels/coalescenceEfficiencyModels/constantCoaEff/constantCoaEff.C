#include "constantCoaEff.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseDynamicMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceEfficiencyModels
{
    defineTypeNameAndDebug(constantCoaEff, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyModel,
        constantCoaEff,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::constantCoaEff::constantCoaEff
(
    const orderedPhasePair& pair,
    const dictionary& dict
)
:
    coalescenceEfficiencyModel(pair, dict),
    K_
    (
        "K",
        dimless,
        dict.subDict("constantEfficiencyCoeffs")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyModels::constantCoaEff::~constantCoaEff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyModels::constantCoaEff::efficiency
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
            K_
        )
    );

    return tE;
}

// ************************************************************************* //

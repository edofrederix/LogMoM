#include "constantCoaEff.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "phaseCompressibleMomentumTransportModel.H"

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
    const dispersedPhaseInterface& pair,
    const dictionary& dict
)
:
    coalescenceEfficiencyModel
    (
        pair,
        dict.subDict(this->type() + this->coeffsDictName_)
    ),
    K_("K", dimless, coeffs_)
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

#include "SugrueLift.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(SugrueLift, 0);
    addToRunTimeSelectionTable(liftModel, SugrueLift, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::SugrueLift::SugrueLift
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedLiftModel(dict, interface),
    aspectRatio_(aspectRatioModel::New(dict.subDict("aspectRatio"), interface))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::SugrueLift::~SugrueLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::SugrueLift::Cl() const
{
    const volScalarField k(continuousTurbulence().k());

    const volScalarField EoH
    (
        interface_.Eo(interface_.dispersed().d()/cbrt(aspectRatio_->E()))
    );

    const volScalarField Wo
    (
        EoH*continuousTurbulence().k()
      / sqr
        (
            max
            (
                mag(interface_.Ur()),
                dimensionedScalar(dimLength/dimTime, 1e-4)
            )
        )
    );

    const volScalarField fWo
    (
        min(0.03, (5.0404 - 5.0781*pow(Wo, 0.0108)))
    );

    const volScalarField fAlpha
    (
        1.0155 - 0.0154*exp(8.0506*interface_.dispersed())
    );

    return fWo*fAlpha;
}

const Foam::phaseCompressible::momentumTransportModel&
Foam::liftModels::SugrueLift::continuousTurbulence() const
{
    return
        interface_.phase1().mesh().lookupObject
        <
            phaseCompressible::momentumTransportModel
        >
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                interface_.continuous().name()
            )
        );
}


// ************************************************************************* //

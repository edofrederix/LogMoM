#include "LogMoMInterfacialModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LogMoMInterfacialModel, 0);

    const Foam::NamedEnum
    <
        LogMoMInterfacialModel::integrationMode,
        3
    >
    Foam::LogMoMInterfacialModel::integrationModeNames_
    {"none", "quadrature", "diameter"};

    const scalarList Foam::LogMoMInterfacialModel::pd_
    (
        IStringStream
        (
            "(0 1)"
        )()
    );

    const scalarList Foam::LogMoMInterfacialModel::qd_
    (
        IStringStream
        (
            "(2)"
        )()
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LogMoMInterfacialModel::LogMoMInterfacialModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    integrationMode_
    (
        integrationModeNames_[word(dict.lookup("mode"))]
    ),
    p_(dict.lookupOrDefault<scalarList>("p", pd_)),
    q_(dict.lookupOrDefault<scalarList>("q", qd_)),
    dispersedInterface_
    (
        interface.modelCast<LogMoMInterfacialModel, dispersedPhaseInterface>()
    )
{
    if (integrationMode_ == integrationMode::quadrature)
    {
        quadrature_.reset
        (
            GaussQuadrature::New
            (
                "GHQ",
                readLabel(dict.lookup("GaussHermite"))
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LogMoMInterfacialModel::~LogMoMInterfacialModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::diameterModels::LogMoM& Foam::LogMoMInterfacialModel::logmom() const
{
    return
        dispersedInterface_.mesh().lookupObject<diameterModels::LogMoM>
        (
            IOobject::groupName
            (
                dispersedInterface_.dispersed().name(),
                diameterModels::LogMoM::typeName
            )
        );
}

Foam::scalar Foam::LogMoMInterfacialModel::gamma() const
{
    return logmom().gamma();
}

// ************************************************************************* //

#include "LogMoMInterfacialModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LogMoMInterfacialModel, 0);

    template<>
    const char*
        NamedEnum<LogMoMInterfacialModel::integrationMode, 3>::names[] =
        {"quadrature", "Riemann", "diameter"};

    const NamedEnum<LogMoMInterfacialModel::integrationMode, 3>
        LogMoMInterfacialModel::integrationModeNames_;

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
    Nr_(dict.lookupOrDefault<label>("Nr", -1)),
    blendVelocity_(dict.lookupOrDefault<Switch>("blendVelocity", false)),
    interface_
    (
        interface.modelCast<LogMoMInterfacialModel, dispersedPhaseInterface>()
    )
{
    if (interface_.dispersed().dPtr()->type() != "threeMomentLogNormal")
    {
        FatalErrorInFunction
            << "This model only works with threeMomentLogNormal as "
            << "diameter model."
            << exit(FatalError);
    }

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

    if (integrationMode_ == integrationMode::Riemann)
    {
        if (Nr_ == -1)
        {
            FatalErrorInFunction
                << "No value set for the number of Riemann bins"
                << endl << abort(FatalError);
        }
    }

    diameterModelPtr_ =
        static_cast<diameterModels::threeMomentLogNormal*>
        (
            &const_cast<autoPtr<diameterModel>&>
            (
                interface_.dispersed().dPtr()
            )()
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LogMoMInterfacialModel::~LogMoMInterfacialModel()
{
    diameterModelPtr_ = nullptr;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::LogMoMInterfacialModel::U
(
    const label gamma
) const
{
    const phaseModel& phase = interface_.dispersed();

    if (gamma == 3)
    {
        return phase.U();
    }
    else
    {
        return phase.mesh().lookupObject<volVectorField>
        (
            IOobject::groupName("U" + Foam::name(gamma), phase.name())
        );
    }
}

Foam::tmp<Foam::volVectorField> Foam::LogMoMInterfacialModel::UBlend
(
    const volScalarField& d,
    const volScalarField& d0,
    const volScalarField& d2,
    const volScalarField& d3,
    const volVectorField& U0,
    const volVectorField& U2,
    const volVectorField& U3
) const
{
    const dimensionedScalar ds(d.dimensions(), 1e-16);

    const volScalarField W0
    (
        (d2-max(d,d0))/(d2-d0+ds)*pos(d2-d)
    );

    const volScalarField W3
    (
        (min(d,d3)-d2)/(d3-d2+ds)*neg0(d2-d)
    );

    const volScalarField W2
    (
        max(1.0-W0,0.0)*pos0(d2-d)
      + max(1.0-W3,0.0)*neg(d2-d)
    );

    return W0*U0 + W2*U2 + W3*U3;
}

// ************************************************************************* //

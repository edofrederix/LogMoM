#include "LogMoM.H"
#include "LogMoMSource.H"
#include "fvm.H"
#include "fvc.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(LogMoM, 0);
    addToRunTimeSelectionTable(diameterModel, LogMoM, dictionary);
}
}

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::diameterModels::LogMoM::limitMoments()
{
    // To assure that the Sauter-mean diameter is in between dMin and dMax, the
    // scaled interfacial area is limited accordingly. To assure realizability
    // (i.e., a real non-negative sigma), we must satisfy the relation lambda >=
    // kappai^3/(36*1e6*pi). This condition can be derived from Eq. (9) in
    // Habiyaremye et al., (2022). This condition is extended to assure that the
    // size distribution width remains between sigmaMin and sigmaMax.

    const scalar pi(constant::mathematical::pi);

    const dimensionedScalar kappaiMin(6.0/dMax_);
    const dimensionedScalar kappaiMax(6.0/dMin_);

    kappai_ = min(max(kappai_, kappaiMin), kappaiMax);
    kappai_.correctBoundaryConditions();

    const volScalarField lambdaMin
    (
        pow(kappai_,3.0)/(36.0*pi*1e6)
      * exp(3.0*sqr(sigmaMin_))
    );

    const volScalarField lambdaMax
    (
        pow(kappai_,3.0)/(36.0*pi*1e6)
      * exp(3.0*sqr(sigmaMax_))
    );

    lambda_ = min(max(lambda_, lambdaMin), lambdaMax);
    lambda_.correctBoundaryConditions();
}

void Foam::diameterModels::LogMoM::correctDistribution()
{
    const scalar pi(constant::mathematical::pi);

    sigma_ =
        sqrt
        (
            log
            (
                max
                (
                    cbrt(36.0*pi*1e6*lambda_)/kappai_,
                    dimensionedScalar(dimless, 1.0)
                )
            )
        );

    dcm_ = 6.0/kappai_*exp(-2.5*sqr(sigma_));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::LogMoM::LogMoM
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(phase.name(), this->typeName),
            phase.time().constant(),
            phase.mesh()
        )
    ),
    diameterModel(diameterProperties, phase),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(dimLength, 0.0)
    ),
    continuousPhasePtr_(nullptr),
    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    kappai_
    (
        IOobject
        (
            IOobject::groupName("kappai", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    dcm_
    (
        IOobject
        (
            IOobject::groupName("dcm", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(dimLength, 0.0)
    ),
    sigma_
    (
        IOobject
        (
            IOobject::groupName("sigma", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(dimless, 0.0)
    ),
    dMin_("dMin", dimLength, diameterProperties),
    dMax_("dMax", dimLength, diameterProperties),
    sigmaMin_
    (
        "sigmaMin",
        dimless,
        diameterProperties.lookupOrDefault<scalar>("sigmaMin", 0)
    ),
    sigmaMax_
    (
        "sigmaMax",
        dimless,
        diameterProperties.lookupOrDefault<scalar>("sigmaMax", 2)
    ),
    sources_(diameterProperties.lookup("sources"), LogMoMSource::iNew(*this)),
    gamma_(3),
    fields_()
{
    // Correct the LogMoM model. Do not limit the moments yet during
    // construction because the boundary update in the limit procedure may
    // require fields that do not yet exist.

    correctDistribution();
    d_ = d(3,2);

    // Add fields to the multivariate convection scheme field table

    fields_.add(kappai_);
    fields_.add(lambda_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::LogMoM::~LogMoM()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseCompressible::momentumTransportModel&
Foam::diameterModels::LogMoM::continuousTurbulence() const
{
    return
        phase().mesh().lookupType<phaseCompressible::momentumTransportModel>
        (
            continuousPhase().name()
        );
}

Foam::dispersedPhaseInterface
Foam::diameterModels::LogMoM::interface() const
{
    return dispersedPhaseInterface(phase(), continuousPhase());
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::LogMoM::d
(
    const scalar p,
    const scalar q
) const
{
    // Calculate the diameter that is related to the moments p and q, see
    // Frederix et al. (2019), Eq. (26). The calculation uses the size
    // distribution parameters because they are correctly limited.

    return min(max(dcm_*exp((p + q)/2.0*sqr(sigma_)), dMin_), dMax_);
}

void Foam::diameterModels::LogMoM::correct()
{
    const volScalarField& alpha = phase();
    const volScalarField& rho = phase().rho();
    const surfaceScalarField& phi = phase().phi();
    const surfaceScalarField& alphaPhi = phase().alphaPhi();

    // Initialise the accumulated source terms to the dilatation effect

    fvScalarMatrix R0(lambda_, inv(dimTime));

    fvScalarMatrix R2
    (
      - fvm::SuSp
        (
            (1.0/3.0)
          * (
                (fvc::ddt(alpha) + fvc::div(phase().alphaPhi()))
              - (fvc::ddt(alpha, rho) + fvc::div(phase().alphaRhoPhi()))/rho
            ),
            kappai_
        )
    );

    // Accumulate the run-time selectable sources

    forAll(sources_, j)
    {
        R0 += sources_[j].R(lambda_);
        R2 += sources_[j].R(kappai_);
    }

    const volScalarField alphar(max(alpha, phase().residualAlpha()));

    const Foam::fvModels& fvModels =
        Foam::fvModels::New(phase().mesh());
    const Foam::fvConstraints& fvConstraints =
        Foam::fvConstraints::New(phase().mesh());

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            phase().mesh(),
            fields_,
            alphaPhi,
            phase().mesh().schemes().div
            (
                "div("+phi.name()+","+alpha.name()+")"
            )
        )
    );

    const dimensionedScalar deltaT(phase().mesh().time().deltaT());

    // Interfacial area concentration equation

    fvScalarMatrix kappaiEqn
    (
        fvm::ddt(alpha, kappai_)
      + mvConvection->fvmDiv(alphaPhi, kappai_)
        ==
        R2
      + fvModels.source(alpha, rho, kappai_)/rho
      - correction
        (
            fvm::Sp
            (
                max(phase().residualAlpha() - alpha, scalar(0))/deltaT,
                kappai_
            )
        )
    );

    kappaiEqn.relax();
    fvConstraints.constrain(kappaiEqn);
    kappaiEqn.solve();
    fvConstraints.constrain(kappai_);

    // Number concentration equation

    fvScalarMatrix lambdaEqn
    (
        fvm::ddt(alpha, lambda_)
      + mvConvection->fvmDiv(alphaPhi, lambda_)
        ==
        R0
      + fvModels.source(alpha, rho, lambda_)/rho
      - correction
        (
            fvm::Sp
            (
                max(phase().residualAlpha() - alpha, scalar(0))/deltaT,
                lambda_
            )
        )
    );

    lambdaEqn.relax();
    fvConstraints.constrain(lambdaEqn);
    lambdaEqn.solve();
    fvConstraints.constrain(lambda_);

    // Limit the moments and then update the distribution from those

    limitMoments();
    correctDistribution();

    // Set the Sauter-mean diameter as the representative size

    d_ = d(3,2);

    // Print some useful numbers

    Info<< type() << ": " << phase().name() << endl << incrIndent;

    Info<< indent << "min/mean/max " << kappai_.name() << " "
        << gMin(kappai_) << " "
        << gAverage(kappai_) << " "
        << gMax(kappai_) << endl;

    Info<< indent << "min/mean/max " << lambda_.name() << " "
        << gMin(lambda_) << " "
        << gAverage(lambda_) << " "
        << gMax(lambda_) << endl;

    Info<< indent << "min/mean/max " << d_.name() << " "
        << gMin(d_) << " "
        << gAverage(d_) << " "
        << gMax(d_) << endl;

    Info<< indent << "min/mean/max " << sigma_.name() << " "
        << gMin(sigma_) << " "
        << gAverage(sigma_) << " "
        << gMax(sigma_) << endl;

    Info<< decrIndent;
}

Foam::scalar Foam::diameterModels::LogMoM::setGamma(const scalar gamma)
{
    const scalar gammaPrev(gamma_);
    gamma_ = gamma;
    return gammaPrev;
}

bool Foam::diameterModels::LogMoM::writeData(Ostream& os) const
{
    return os.good();
}

bool Foam::diameterModels::LogMoM::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    diameterProperties().lookup("dMin") >> dMin_;
    diameterProperties().lookup("dMax") >> dMax_;

    sigmaMin_ = diameterProperties().lookupOrDefault<scalar>("sigmaMin", 0);
    sigmaMax_ = diameterProperties().lookupOrDefault<scalar>("sigmaMax", 2);

    PtrList<LogMoMSource>
    (
        diameterProperties().lookup("sources"),
        LogMoMSource::iNew(*this)
    ).transfer(sources_);

    return true;
}

// ************************************************************************* //

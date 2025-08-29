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

void Foam::diameterModels::LogMoM::correctLimitedMoments()
{
    // To assure that the Sauter-mean diameter is in between dMin and dMax, the
    // scaled interfacial area is limited accordingly. To assure realizability
    // (i.e., a real non-negative sigma), we must satisfy the relation lambda >=
    // kappa^3/(36*1e6*pi). This condition can be derived from Eq. (9) in
    // Habiyaremye et al., (2022). This condition is extended to assure that the
    // size distribution width remains between sigmaMin and sigmaMax.

    const volScalarField& alpha = phase();
    const dimensionedScalar residualAlpha(phase().residualAlpha());

    const scalar pi(constant::mathematical::pi);

    const volScalarField alphar(max(alpha,residualAlpha));

    const dimensionedScalar kappaMin(6.0/dMax_);
    const dimensionedScalar kappaMax(6.0/dMin_);

    kappa_ = min(max(A_/alphar, kappaMin), kappaMax);
    kappa_.correctBoundaryConditions();

    const volScalarField lambdaMin
    (
        pow(kappa_,3.0)/(36.0*pi*1e6)
      * exp(3.0*sqr(sigmaMin_))
    );

    const volScalarField lambdaMax
    (
        pow(kappa_,3.0)/(36.0*pi*1e6)
      * exp(3.0*sqr(sigmaMax_))
    );

    lambda_ = min(max(N_/alphar, lambdaMin), lambdaMax);
    lambda_.correctBoundaryConditions();
}

void Foam::diameterModels::LogMoM::correctDistribution()
{
    // Update the width and count median diameter of the size distribution from
    // the limited moments

    const scalar pi(constant::mathematical::pi);

    sigma_ =
        sqrt
        (
            log
            (
                max
                (
                    cbrt(36.0*pi*1e6*lambda_)/kappa_,
                    dimensionedScalar(dimless, 1.0)
                )
            )
        );

    dcm_ = 6.0/kappa_*exp(-2.5*sqr(sigma_));
}

void Foam::diameterModels::LogMoM::correctRepresentativeDiameter()
{
    d_ = d(p_, q_);
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
    p_(diameterProperties.lookupOrDefault<scalar>("p", 3)),
    q_(diameterProperties.lookupOrDefault<scalar>("q", 2)),
    continuousPhasePtr_(nullptr),
    N_
    (
        IOobject
        (
            IOobject::groupName("N", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    A_
    (
        IOobject
        (
            IOobject::groupName("A", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(N_.dimensions(), 0.0)
    ),
    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            phase.mesh().time().name(),
            phase.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(A_.dimensions(), 0.0)
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
    nCorr_(diameterProperties.lookupOrDefault<label>("nCorr", 2)),
    sources_(diameterProperties.lookup("sources"), LogMoMSource::iNew(*this)),
    gamma_(3),
    fields_()
{
    // Correct the LogMoM model

    correctLimitedMoments();
    correctDistribution();
    correctRepresentativeDiameter();

    // Add fields to the multivariate convection scheme field table

    fields_.add(kappa_);
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
    // Calculate the diameter that is related to moment p and q, see Frederix et
    // al. (2019), Eq. (26).

    return min(max(dcm_*exp((p + q)/2.0*sqr(sigma_)), dMin_), dMax_);
}

void Foam::diameterModels::LogMoM::correct()
{
    const volScalarField& alpha = phase();
    const volScalarField& rho = phase().rho();
    const surfaceScalarField& phi = phase().phi();
    const surfaceScalarField& alphaPhi = phase().alphaPhi();

    // Initialise the accumulated source terms to the dilatation effect

    fvScalarMatrix R0(N_, inv(dimTime));

    fvScalarMatrix R2
    (
      - fvm::SuSp
        (
            (1.0/3.0)
          * (
                (fvc::ddt(alpha) + fvc::div(phase().alphaPhi()))
              - (fvc::ddt(alpha, rho) + fvc::div(phase().alphaRhoPhi()))/rho
            ),
            A_
        )
    );

    // Accumulate the run-time selectable sources

    forAll(sources_, j)
    {
        R0 += sources_[j].R(N_);
        R2 += sources_[j].R(A_);
    }

    const volScalarField alphar(max(alpha, phase().residualAlpha()));

    const Foam::fvModels& fvModels =
        Foam::fvModels::New(phase().mesh());
    const Foam::fvConstraints& fvConstraints =
        Foam::fvConstraints::New(phase().mesh());

    for (int corr = 0; corr < nCorr_; corr++)
    {
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

        fvScalarMatrix AEqn
        (
            fvm::ddt(A_)
          + mvConvection->fvcDiv(alphaPhi, kappa_)
            ==
            R2
          + fvModels.source(alpha, rho, A_)/(alphar*rho)
        );

        AEqn.relax();
        fvConstraints.constrain(AEqn);

        A_ = AEqn.H()/AEqn.A();
        A_.correctBoundaryConditions();

        fvConstraints.constrain(A_);

        fvScalarMatrix NEqn
        (
            fvm::ddt(N_)
          + mvConvection->fvcDiv(alphaPhi, lambda_)
            ==
            R0
          + fvModels.source(alpha, rho, N_)/(alphar*rho)
        );

        NEqn.relax();
        fvConstraints.constrain(NEqn);

        N_ = NEqn.H()/NEqn.A();
        N_.correctBoundaryConditions();

        fvConstraints.constrain(N_);

        // Correct limited moments and update the main moments from those

        correctLimitedMoments();

        A_ = kappa_*alpha;
        N_ = lambda_*alpha;

        A_.correctBoundaryConditions();
        N_.correctBoundaryConditions();
    }

    correctDistribution();
    correctRepresentativeDiameter();

    Info<< A_.name() << ", min, max = "
        << gAverage(A_) << " "
        << gMin(A_) << " "
        << gMax(A_) << endl;

    Info<< N_.name() << ", min, max = "
        << gAverage(N_) << " "
        << gMin(N_) << " "
        << gMax(N_) << endl;
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

    p_ = diameterProperties().lookupOrDefault<scalar>("p", 3);
    q_ = diameterProperties().lookupOrDefault<scalar>("q", 2);

    diameterProperties().lookup("dMin") >> dMin_;
    diameterProperties().lookup("dMax") >> dMax_;

    sigmaMin_ = diameterProperties().lookupOrDefault<scalar>("sigmaMin", 0);
    sigmaMax_ = diameterProperties().lookupOrDefault<scalar>("sigmaMax", 2);

    nCorr_ = diameterProperties().lookupOrDefault<label>("nCorr", 2);

    PtrList<LogMoMSource>
    (
        diameterProperties().lookup("sources"),
        LogMoMSource::iNew(*this)
    ).transfer(sources_);

    return true;
}

// ************************************************************************* //

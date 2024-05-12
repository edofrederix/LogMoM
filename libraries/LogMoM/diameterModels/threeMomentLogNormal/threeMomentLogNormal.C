#include "threeMomentLogNormal.H"
#include "phaseSystem.H"
#include "fvm.H"
#include "fvc.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "GaussQuadrature.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(threeMomentLogNormal, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        threeMomentLogNormal,
        dictionary
    );
}
}

// * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * *//

void Foam::diameterModels::threeMomentLogNormal::correctLimitedScaledMoments()
{
    // In the limit of small alpha, we want to move to a mono-dispersed size
    // distribution at size dMin. This is done by the blending function F. To
    // assure that the Sauter mean diameter is in between dMin and dMax, the
    // scaled interfacial area is limited accordingly. Furthermore, to assure
    // realizability (i.e., a real non-negative sigma), we must satisfy the
    // relation lambda >= kappa^3/(36*1e6*pi). This condition can be derived
    // from Eq. (9) in Habiyaremye et al., (2022).

    const volScalarField& alpha = phase();
    const dimensionedScalar residualAlpha(phase().residualAlpha());

    const volScalarField F
    (
        max(residualAlpha - alpha, scalar(0))/residualAlpha
    );

    const scalar pi(constant::mathematical::pi);

    const volScalarField alphar(max(alpha,residualAlpha));

    kappa_ =
        min
        (
            max
            (
                6.0/dMin_*F + A_/alphar*(1.0-F),
                6.0/dMax_
            ),
            6.0/dMin_
        );

    lambda_ =
        max
        (
            N_/alphar,
            pow(kappa_,3.0)/(36.0*pi*1e6)
        );
}

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::threeMomentLogNormal::dsm() const
{
    return max(6/max(kappa_, 6/dMax_), dMin_);
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::threeMomentLogNormal::d
(
    const scalar p,
    const scalar q
) const
{
    const scalar pi(constant::mathematical::pi);

    tmp<volScalarField> td
    (
        new volScalarField
        (
            IOobject
            (
                "d(" + Foam::name(p) + "," + Foam::name(q) + ")",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimLength, 0)
        )
    );

    volScalarField& d = td.ref();

    d = min
    (
        max
        (
            pow(6.0,(p+q-2)/3)
          * pow(pi*lambda_*1e6,(p+q-5)/6)
          / pow(kappa_,(p+q-3)/2),
            dMin_
        ),
        dMax_
    );

    return td;
}

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::threeMomentLogNormal::sigma() const
{
    const scalar pi(constant::mathematical::pi);

    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField& s = tsigma.ref();

    s =
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

    return tsigma;
}


void Foam::diameterModels::threeMomentLogNormal::updateCoalescenceSources
(
    const volScalarField& dcm,
    const volScalarField& s
)
{
    const scalar pi(constant::mathematical::pi);

    const volScalarField alpha(phase());

    // Gauss-Hermite quadratures for integration of coalescence sources

    autoPtr<GaussQuadrature> GHQPtr
    (
        GaussQuadrature::New
        (
            "GHQ",
            readLabel
            (
                diameterProperties().subDict("coalescence")
               .lookup("GaussHermite")
            )
        )
    );

    const GaussQuadrature& GHQ = GHQPtr();

    // Diameter fields corresponding to quadrature nodes

    PtrList<volScalarField> d(GHQ.size());

    forAll(d, i)
    {
        const scalar x(GHQ.x()[i]);

        d.set
        (
            i,
            new volScalarField
            (
                exp(x*sqrt(2.0)*s)*dcm
            )
        );
    }

    volScalarField NCoa
    (
        IOobject
        (
             "NCoa",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(dimVolume/dimTime, Zero)
    );

    volScalarField ACoa
    (
        IOobject
        (
             "ACoa",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(pow5(dimLength)/dimTime, Zero)
    );

    for (label i = 0; i < GHQ.size(); i++)
    {
        const scalar wi(GHQ.w()[i]);

        for (label j = 0; j <= i; j++)
        {
            const scalar wj(GHQ.w()[j]);

            const volScalarField K
            (
                coaEffModel_->efficiency(d[i],d[j])
              * coaFreqModel_->frequency(d[i],d[j])
            );

            NCoa += -K*wi*wj*(i == j ? 1.0 : 2.0);

            ACoa +=
                (
                    cbrt(sqr(pow3(d[i]) + pow3(d[j])))
                  - sqr(d[i])
                  - sqr(d[j])
                )
              * K*wi*wj*pi
              * (i == j ? 1.0 : 2.0);
        }
    }

    NCoaRate_ = -0.5/pi*alpha*lambda_*1e6*NCoa;
    ACoaRate_ = -0.5/pi*alpha*sqr(1e6*lambda_)/kappa_*ACoa;
}

void Foam::diameterModels::threeMomentLogNormal::updateBreakupSources
(
    const volScalarField& dcm,
    const volScalarField& s
)
{
    const scalar pi(constant::mathematical::pi);

    const volScalarField alpha(phase());

    // Gauss-Hermite quadratures for integration of breakup sources

    autoPtr<GaussQuadrature> GHQPtr
    (
        GaussQuadrature::New
        (
            "GHQ",
            readLabel
            (
                diameterProperties().subDict("breakup")
               .lookup("GaussHermite")
            )
        )
    );

    // Gauss-Legendre quadratures for integration of breakup sources

    autoPtr<GaussQuadrature> GLQPtr
    (
        GaussQuadrature::New
        (
            "GLQ",
            readLabel
            (
                diameterProperties().subDict("breakup")
               .lookup("GaussLegendre")
            )
        )
    );

    const GaussQuadrature& GHQ = GHQPtr();
    const GaussQuadrature& GLQ = GLQPtr();

    // Diameter field corresponding to quadrature nodes

    PtrList<volScalarField> d1(GHQ.size());

    // Partial break-up rates

    PtrList<volScalarField> Blist(ceil(GLQ.size()/2.0));

    forAll(d1, i)
    {
        const scalar x(GHQ.x()[i]);

        d1.set
        (
            i,
            new volScalarField
            (
                exp(x*sqrt(2.0)*s)*dcm
            )
        );
    }

    volScalarField NBreak
    (
        IOobject
        (
             "NBreak",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), Zero)
    );

    volScalarField ABreak
    (
        IOobject
        (
             "ABreak",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(sqr(dimLength)/dimTime, Zero)
    );

    volScalarField B
    (
        IOobject
        (
             "B",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(inv(dimVolume*dimTime), Zero)
    );

    for (label i = 0; i < GHQ.size(); i++)
    {
        const scalar wi(GHQ.w()[i]);

        for (label j = 0; j < GLQ.size(); j++)
        {
            const scalar x(GLQ.x()[j]);

            // Integration variable for Gauss-Legendre quadrature

            const volScalarField w((x+1)*pow3(d1[i])/2.0);
            const volScalarField d2(cbrt(w));

            if (x <= 0)
            {
                Blist.set
                (
                    j,
                    new volScalarField
                    (
                        breakupModel_->binaryRate(d1[i], d2)
                    )
                );
                B = Blist[j];
            }
            else
            {
                B = Blist[GLQ.size()-1-j];
            }

            const scalar wj(GLQ.w()[j]);

            NBreak +=
                B*(pi/2.0)*pow3(d1[i])
              * 0.5*wj*wi*(1.0/3.0);

            ABreak +=
                B*(pi/2.0)*pow3(d1[i])
              * (sqr(d2)-sqr(d1[i])/2.0)
              * wj*wi*(1.0/3.0);
        }
    }

    NBreakRate_ = 1.0/(2.0*sqrt(pi))*NBreak;
    ABreakRate_ = sqrt(pi)/2.0*1e6*lambda_/kappa_*ABreak;
}

void Foam::diameterModels::threeMomentLogNormal::readModels()
{
    if (coalescence_ && coaEffModel_.empty() && coaFreqModel_.empty())
    {
        coaEffModel_ =
            coalescenceEfficiencyModel::New
            (
                interface(),
                diameterProperties().subDict("coalescence")
            );

        coaFreqModel_ =
            coalescenceFrequencyModel::New
            (
                interface(),
                diameterProperties().subDict("coalescence")
            );
    }

    if (breakup_ && breakupModel_.empty())
    {
        breakupModel_ =
            breakupModel::New
            (
                interface(),
                diameterProperties().subDict("breakup")
            );
    }
}

void Foam::diameterModels::threeMomentLogNormal::updateSources()
{
    if (coalescence_ || breakup_)
    {
        readModels();

        const volScalarField s(sigma());
        const volScalarField dcm(d(0.0,0.0));

        if (coalescence_)
        {
            updateCoalescenceSources(dcm, s);
        }

        if (breakup_)
        {
            updateBreakupSources(dcm, s);
        }
    }

    const phaseModel& phase = this->phase();

    // Break-up and coalescence

    NSource_ = NBreakRate_ - NCoaRate_;
    ASource_ = ABreakRate_ - ACoaRate_;

    // Dilation

    ASource_ -=
        (1.0/3.0)
      * (
            (
                fvc::ddt(phase)
              + fvc::div(phase.alphaPhi())
            )
          - (
                fvc::ddt(phase, phase.rho()())
              + fvc::div(phase.alphaRhoPhi())
            )
          / phase.rho()
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::threeMomentLogNormal::threeMomentLogNormal
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    diameterModel(dict, phase),
    mesh_(phase.mesh()),
    p_(readScalar(dict.lookup("p"))),
    q_(readScalar(dict.lookup("q"))),
    N_
    (
        IOobject
        (
            IOobject::groupName("N", phase.name()),
            mesh_.time().timeName(),
            mesh_,
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
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(N_.dimensions(), 0.0)
    ),
    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(A_.dimensions(), 0.0)
    ),
    dMax_("dMax", dimLength, dict),
    dMin_("dMin", dimLength, dict),
    nCorr_(dict.lookupOrDefault<label>("nCorr", 2)),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, 0.0)
    ),
    sigma_
    (
        IOobject
        (
            IOobject::groupName("sigma", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    ),
    coalescence_(dict.subDict("coalescence").lookup("active")),
    breakup_(dict.subDict("breakup").lookup("active")),
    NCoaRate_
    (
        IOobject
        (
            IOobject::groupName("NCoaRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    ACoaRate_
    (
        IOobject
        (
            IOobject::groupName("ACoaRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    NBreakRate_
    (
        IOobject
        (
            IOobject::groupName("NBreakRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    ABreakRate_
    (
        IOobject
        (
            IOobject::groupName("ABreakRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    NSource_
    (
        IOobject
        (
            IOobject::groupName("NSource", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(inv(dimTime), 0)
    ),
    ASource_
    (
        IOobject
        (
            IOobject::groupName("ASource", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phase.mesh(),
        dimensionedScalar(inv(dimTime), 0)
    )
{
    correctLimitedScaledMoments();

    d_ = dsm();
    sigma_ = sigma();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::threeMomentLogNormal::~threeMomentLogNormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::threeMomentLogNormal::correct()
{
    const phaseModel& phase = this->phase();
    const volScalarField& alpha = phase;
    const surfaceScalarField& alphaPhi = phase.alphaPhi();

    // Try to find poly-celerity fluxes. Otherwise reduce to mono-celerity ones.

    const word alphaPhi0Name(IOobject::groupName("alphaPhi0", phase.name()));
    const surfaceScalarField& alphaPhi0
    (
        mesh_.foundObject<surfaceScalarField>(alphaPhi0Name)
      ? mesh_.lookupObjectRef<surfaceScalarField>(alphaPhi0Name)
      : alphaPhi
    );

    const word alphaPhi2Name(IOobject::groupName("alphaPhi2", phase.name()));
    const surfaceScalarField& alphaPhi2
    (
        mesh_.foundObject<surfaceScalarField>(alphaPhi2Name)
      ? mesh_.lookupObjectRef<surfaceScalarField>(alphaPhi2Name)
      : alphaPhi
    );

    updateSources();

    const Foam::fvModels& fvModels(Foam::fvModels::New(phase.mesh()));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(phase.mesh())
    );

    for (int corr = 0; corr < nCorr_; corr++)
    {
        // Use non-limited scaled moments, for stability

        kappa_ = A_/max(alpha, phase.residualAlpha());
        lambda_ = N_/max(alpha, phase.residualAlpha());

        fvScalarMatrix AEqn
        (
            fvm::ddt(A_)
          + fvc::div(alphaPhi2, kappa_)
          ==
          - fvm::SuSp(-ASource_, A_)
          + fvModels.source(A_)
        );

        AEqn.relax();
        fvConstraints.constrain(AEqn);

        A_ = AEqn.H()/AEqn.A();
        A_.correctBoundaryConditions();

        fvScalarMatrix NEqn
        (
            fvm::ddt(N_)
          + fvc::div(alphaPhi0, lambda_)
          ==
          - fvm::SuSp(-NSource_, N_)
          + fvModels.source(N_)
        );

        NEqn.relax();
        fvConstraints.constrain(NEqn);

        N_ = NEqn.H()/NEqn.A();
        N_.correctBoundaryConditions();
    }

    correctLimitedScaledMoments();

    d_ = dsm();
    sigma_ = sigma();
}

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::threeMomentLogNormal::d() const
{
    return d_;
}

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::threeMomentLogNormal::a() const
{
    // Instead of returning the interfacial area directly, we return the scaled
    // one which is limited

    return phase()*kappa_;
}

bool Foam::diameterModels::threeMomentLogNormal::read
(
    const dictionary& phaseProperties
)
{
    diameterModel::read(phaseProperties);

    diameterProperties().lookup("dMax") >> dMax_;
    diameterProperties().lookup("dMin") >> dMin_;

    return true;
}

const Foam::phaseModel&
Foam::diameterModels::threeMomentLogNormal::continuousPhase() const
{
    return
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName
            (
                "alpha",
                diameterProperties().lookup("continuousPhase")
            )
        );
}

Foam::dispersedPhaseInterface
Foam::diameterModels::threeMomentLogNormal::interface() const
{
    return dispersedPhaseInterface(phase(), continuousPhase());
}

const Foam::phaseCompressible::momentumTransportModel&
Foam::diameterModels::threeMomentLogNormal::continuousTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressible::momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                continuousPhase().name()
            )
        );
}

const Foam::phaseCompressible::momentumTransportModel&
Foam::diameterModels::threeMomentLogNormal::dispersedTurbulence() const
{
    return
        mesh_.lookupObject<phaseCompressible::momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                phase().name()
            )
        );
}

// ************************************************************************* //

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

    template<>
    const char*
        NamedEnum<threeMomentLogNormal::closingMomentType, 2>::names[] =
        {"interfacialArea", "squaredVolume"};

    const NamedEnum<threeMomentLogNormal::closingMomentType, 2>
        threeMomentLogNormal::closingMomentTypeNames_;
}
}

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * *//

void Foam::diameterModels::threeMomentLogNormal::limitMoments(bool corrBou)
{
    const scalar pi(constant::mathematical::pi);

    // Limit lambda based on the diameter of average mass

    lambda_ =
        min
        (
            max(lambda_, 6.0/(pi*1e6*pow(dMax_,3))),
            6.0/(pi*1e6*pow(dMin_,3))
        );

    if (corrBou)
    {
        lambda_.correctBoundaryConditions();
    }

    // Limit kappa based on the Sauter mean diameter

    kappa_ =
        min
        (
            max(kappa_, 6.0/dMax_),
            6.0/dMin_
        );

    if (corrBou)
    {
        kappa_.correctBoundaryConditions();
    }

    // Limit beta based on d_36

    beta_ =
        min
        (
            max(beta_, 6.0*pow(dMin_,3)/pi),
            6.0*pow(dMax_,3)/pi
        );

    if (corrBou)
    {
        beta_.correctBoundaryConditions();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::threeMomentLogNormal::dsm() const
{
    return d(p_,q_);
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

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        // Compute the diameter from the interfacial area concentration

        d = pow(6.0,(p+q-2)/3)
          * pow(pi*lambda_*1e6,(p+q-5)/6)
          / pow(kappa_,(p+q-3)/2);
    }
    else
    {
        // Compute the diameter from the squared volume concentration

        d = pow(pi/6.0,(p+q-6)/9)
          * pow(beta_,(p+q-3)/18)
          * pow(1e6*lambda_,(p+q-9)/18);
    }

    return min(max(d,dMin_),dMax_);
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

    const dimensionedScalar unityDimless(dimless,1.0);

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        s =
            sqrt
            (
                log
                (
                    max
                    (
                        cbrt(36.0*pi*1e6*lambda_)/kappa_,
                        unityDimless
                    )
                )
            );
    }
    else
    {
        s =
            sqrt
            (
                1.0/9.0
              * log
                (
                    max
                    (
                        beta_*1e6*lambda_*sqr(pi/6.0),
                        unityDimless
                    )
                )
            );
    }

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

    volScalarField lambdaCoa
    (
        IOobject
        (
             "lambdaCoa",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(dimVolume/dimTime, Zero)
    );

    volScalarField kappaCoa
    (
        IOobject
        (
             "kappaCoa",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(pow5(dimLength)/dimTime, Zero)
    );

    volScalarField betaCoa
    (
        IOobject
        (
             "betaCoa",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(pow(dimLength,9)/dimTime, Zero)
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

            lambdaCoa += -K*wi*wj
                       * (i == j ? 1.0 : 2.0);

            if (closingMoment_ == closingMomentType::interfacialArea)
            {
                kappaCoa +=
                    (
                        cbrt(sqr(pow3(d[i]) + pow3(d[j])))
                      - sqr(d[i])
                      - sqr(d[j])
                    )
                  * K*wi*wj*pi
                  * (i == j ? 1.0 : 2.0);
            }
            else
            {
                betaCoa +=
                    2.0*pow3(d[i])*pow3(d[j])
                  * K*wi*wj
                  * (i == j ? 1.0 : 2.0);
            }
        }
    }

    lambdaCoaRate_ = -0.5/pi*sqr(alpha)*1e6*lambda_*lambdaCoa;

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        kappaCoaRate_ =
          - 0.5/pi
          * sqr(alpha*1e6*lambda_)
          * kappaCoa
          / kappa_;
    }
    else
    {
        betaCoaRate_ =
          - 0.5/pi
          * sqr(alpha*1e6*lambda_)
          * betaCoa
          / beta_;
    }
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

    volScalarField lambdaBreak
    (
        IOobject
        (
             "lambdaBreak",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), Zero)
    );

    volScalarField kappaBreak
    (
        IOobject
        (
             "kappaBreak",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(sqr(dimLength)/dimTime, Zero)
    );

    volScalarField betaBreak
    (
        IOobject
        (
             "betaBreak",
             mesh_.time().timeName(),
             mesh_
        ),
        mesh_,
        dimensionedScalar(pow6(dimLength)/dimTime, Zero)
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

            lambdaBreak +=
                B*(pi/2.0)*pow3(d1[i])
              * 0.5*wj*wi*(1.0/3.0);

            if (closingMoment_ == closingMomentType::interfacialArea)
            {
                kappaBreak +=
                    B*(pi/2.0)*pow3(d1[i])
                  * (sqr(d2)-sqr(d1[i])/2.0)
                  * wj*wi*(1.0/3.0);
            }
            else
            {
                betaBreak +=
                    B*(pi/2.0)*pow3(d1[i])
                  * (pow6(d2)-pow6(d1[i])/2.0)
                  * wj*wi*(1.0/3.0);
            }
        }
    }

    lambdaBreakRate_ = alpha/(2.0*sqrt(pi))*lambdaBreak;

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        kappaBreakRate_ =
            alpha*sqrt(pi)/2.0*1e6*lambda_
          * max(kappaBreak, dimensionedScalar(sqr(dimLength)/dimTime, Zero))
          / kappa_;
    }
    else
    {
        betaBreakRate_ =
            alpha/(2.0*sqrt(pi))*1e6*lambda_*betaBreak
          / beta_;
    }
}

void Foam::diameterModels::threeMomentLogNormal::readModels()
{
    const phaseModel& continuousPhase =
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName
            (
                "alpha",
                diameterProperties().lookup("continuousPhase")
            )
        );

    const dispersedPhaseInterface pair(phase(), continuousPhase);

    if (coalescence_ && coaEffModel_.empty() && coaFreqModel_.empty())
    {
        coaEffModel_ =
            coalescenceEfficiencyModel::New
            (
                pair,
                diameterProperties().subDict("coalescence")
            );

        coaFreqModel_ =
            coalescenceFrequencyModel::New
            (
                pair,
                diameterProperties().subDict("coalescence")
            );
    }

    if (breakup_ && breakupModel_.empty())
    {
        breakupModel_ =
            breakupModel::New
            (
                pair,
                diameterProperties().subDict("breakup")
            );
    }
}

void Foam::diameterModels::threeMomentLogNormal::updateModels()
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
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::threeMomentLogNormal::threeMomentLogNormal
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    mesh_(phase.mesh()),
    p_(readScalar(diameterProperties.lookup("p"))),
    q_(readScalar(diameterProperties.lookup("q"))),
    closingMoment_
    (
        closingMomentTypeNames_[word(diameterProperties.lookup("closingMoment"))]
    ),
    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    beta_
    (
        IOobject
        (
            IOobject::groupName("beta", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    dMax_("dMax", dimLength, diameterProperties),
    dMin_("dMin", dimLength, diameterProperties),
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
        dsm()
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
        sigma()
    ),
    coalescence_(diameterProperties.subDict("coalescence").lookup("active")),
    breakup_(diameterProperties.subDict("breakup").lookup("active")),
    lambdaCoaRate_
    (
        IOobject
        (
            IOobject::groupName("lambdaCoaRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    kappaCoaRate_
    (
        IOobject
        (
            IOobject::groupName("kappaCoaRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    betaCoaRate_
    (
        IOobject
        (
            IOobject::groupName("betaCoaRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    lambdaBreakRate_
    (
        IOobject
        (
            IOobject::groupName("lambdaBreakRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    kappaBreakRate_
    (
        IOobject
        (
            IOobject::groupName("kappaBreakRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    ),
    betaBreakRate_
    (
        IOobject
        (
            IOobject::groupName("betaBreakRate", phase.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(inv(dimTime), 0)
    )
{
    limitMoments(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::threeMomentLogNormal::~threeMomentLogNormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::threeMomentLogNormal::correct()
{
    const scalar pi(constant::mathematical::pi);

    const phaseModel& phase = this->phase();
    const surfaceScalarField& alphaPhi = phase.alphaPhi();
    const dimensionedScalar& residualAlpha = phase.residualAlpha();

    kappa_.correctBoundaryConditions();
    lambda_.correctBoundaryConditions();

    volScalarField R
    (
        IOobject
        (
            "divU",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, 0)
    );

    if (phase.divU().valid())
    {
        R = -phase.divU();
    }

    updateModels();

    const Foam::fvModels& fvModels(Foam::fvModels::New(phase.mesh()));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(phase.mesh())
    );

    // Note: the N-equation is formulated in terms of lambda (lambda =
    // N/alpha/1e6), however, it remains proportional to ddt(N) and not
    // ddt(lambda). The factor 1e6 drops out everywhere.

    fvScalarMatrix NEqn
    (
        fvm::ddt(phase, lambda_)
      + fvm::div(alphaPhi, lambda_)
      ==
      - fvm::SuSp(lambdaCoaRate_, lambda_)
      - fvm::SuSp(-lambdaBreakRate_, lambda_)
      + residualAlpha
      * (
            fvc::ddt(lambda_)
          - fvm::ddt(lambda_)
        )
      + fvModels.source(phase, lambda_)
    );

    if (closingMoment_ == closingMomentType::interfacialArea)
    {

        // Note: the ai-equation is formulated in terms of kappa (kappa =
        // ai/alpha), however, it remains proportional to ddt(ai) and not
        // ddt(kappa)

        fvScalarMatrix aiEqn
        (
            fvm::ddt(phase, kappa_)
          + fvm::div(alphaPhi, kappa_)
          ==
          - fvm::SuSp(kappaCoaRate_, kappa_)
          - fvm::SuSp(-kappaBreakRate_, kappa_)
          - fvm::SuSp(2.0/3.0*R, kappa_)
          + residualAlpha
          * (
                fvc::ddt(kappa_)
              - fvm::ddt(kappa_)
            )
          + fvModels.source(phase, kappa_)
        );

        aiEqn.relax();

        fvConstraints.constrain(aiEqn);

        aiEqn.solve();
    }
    else
    {

        // Note: the B-equation is formulated in terms of beta (beta = B/alpha),
        // however, it remains proportional to ddt(B) and not ddt(beta)

        fvScalarMatrix BEqn
        (
            fvm::ddt(phase, beta_)
          + fvm::div(alphaPhi, beta_)
          ==
          - fvm::SuSp(betaCoaRate_, beta_)
          - fvm::SuSp(-betaBreakRate_, beta_)
          - fvm::SuSp(2.0*R, beta_)
          + residualAlpha
          * (
                fvc::ddt(beta_)
              - fvm::ddt(beta_)
            )
          + fvModels.source(phase, beta_)
        );

        BEqn.relax();

        fvConstraints.constrain(BEqn);

        BEqn.solve();
    }

    // Solve the N-equation after the ai- or B-equation has been formulated and
    // solved, for consistency

    NEqn.relax();

    fvConstraints.constrain(NEqn);

    NEqn.solve();

    limitMoments();

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        beta_ = pi*pow(6,8)*sqr(1e6*lambda_)/pow(kappa_, 9);
    }
    else
    {
        kappa_ =
            pow
            (
                pi*pow(6,8)*sqr(1e6*lambda_)/beta_,
                1.0/9.0
            );
    }

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

// ************************************************************************* //

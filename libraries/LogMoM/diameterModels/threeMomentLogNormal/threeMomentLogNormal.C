#include "threeMomentLogNormal.H"
#include "phaseSystem.H"
#include "fvm.H"
#include "fvc.H"
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

    const dimensionedScalar lambdaSmall(lambda_.dimensions(), SMALL);
    const dimensionedScalar kappaSmall(kappa_.dimensions(), pow(SMALL,3));
    const dimensionedScalar betaSmall(beta_.dimensions(), pow(SMALL,7));

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        // Compute the diameter from the interfacial area concentration

        d = max
            (
                min
                (
                    pow(6.0, (p+q-2)/3)
                  * pow(pi*max(1e6*lambda_,lambdaSmall), (p+q-5)/6)
                  / pow(max(kappa_,kappaSmall), (p+q-3)/2),
                    dMax_
                ),
                dMin_
            );
    }
    else
    {
        // Compute the diameter from the squared volume concentration

        const dimensionedScalar betaSmall(beta_.dimensions(), pow(SMALL,7));

        d = max
            (
                min
                (
                    pow(pi/6.0, (p+q-6)/9)
                  * pow(max(beta_,betaSmall), (p+q-3)/18)
                  * pow(max(1e6*lambda_,lambdaSmall), (p+q-9)/18),
                    dMax_
                ),
                dMin_
            );
    }

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

    const dimensionedScalar unityDimless(dimless,1.0);

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        const dimensionedScalar smallKappa(kappa_.dimensions(), pow(SMALL,3));

        s =
            sqrt
            (
                log
                (
                    max
                    (
                        cbrt(36.0*pi*1e6*lambda_)/max(kappa_,smallKappa),
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
                    max(beta_*1e6*lambda_*sqr(pi/6.0), unityDimless)
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
          / max
            (
                kappa_,
                dimensionedScalar
                (
                     inv(dimLength),
                     pow(SMALL,3)
                )
            );
    }
    else
    {
        betaCoaRate_ =
          - 0.5/pi
          * sqr(alpha*1e6*lambda_)
          * betaCoa
          / max
            (
                beta_,
                dimensionedScalar
                (
                     dimVolume,
                     pow(SMALL,7)
                )
            );
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
          / max
            (
                kappa_,
                dimensionedScalar
                (
                    inv(dimLength),
                    pow(SMALL,3)
                )
            );
    }
    else
    {
        betaBreakRate_ =
            alpha/(2.0*sqrt(pi))*1e6*lambda_*betaBreak
          / max
            (
                beta_,
                dimensionedScalar
                (
                    dimVolume,
                    pow(SMALL,7)
                )
            );
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

    const orderedPhasePair pair(phase(), continuousPhase);

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
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::threeMomentLogNormal::~threeMomentLogNormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::threeMomentLogNormal::correct()
{
    const scalar pi(constant::mathematical::pi);

    const phaseModel& phase = this->phase();

    const phaseModel& continuousPhase =
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName
            (
                "alpha",
                diameterProperties().lookup("continuousPhase")
            )
        );

    const surfaceScalarField& alphaPhi = phase.alphaPhi();

    const volVectorField U(phase.U());
    const volVectorField Uc(continuousPhase.U());

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
        );

        aiEqn.relax();
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
        );

        BEqn.relax();
        BEqn.solve();
    }

    // Solve the N-equation after the ai- or B-equation has been formulated and
    // solved, for consistency

    NEqn.relax();
    NEqn.solve();

    if (closingMoment_ == closingMomentType::interfacialArea)
    {
        const dimensionedScalar kappaSmall(kappa_.dimensions(), pow(SMALL,3));

        beta_ =
            pi*pow(6,8)*sqr(1e6*lambda_)
          / pow(max(kappa_, kappaSmall), 9);
    }
    else
    {
        const dimensionedScalar betaSmall(beta_.dimensions(), pow(SMALL,7));

        kappa_ =
            pow
            (
                pi*pow(6,8)*sqr(1e6*lambda_)
              / max(beta_, betaSmall),
                1.0/9.0
            );
    }


    lambda_.max(0.0);
    kappa_.max(0.0);
    beta_.max(0.0);
    beta_.min(1e9);

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

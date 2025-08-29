#include "LogMoMCoalescence.H"
#include "LogMoM.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace LogMoMSources
{
    defineTypeNameAndDebug(LogMoMCoalescence, 0);
    addToRunTimeSelectionTable(LogMoMSource, LogMoMCoalescence, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::LogMoMSources::LogMoMCoalescence::LogMoMCoalescence
(
    const LogMoM& logmom,
    const dictionary& dict
)
:
    LogMoMSource(logmom),
    GHQ_
    (
        GaussQuadrature::New
        (
            "GHQ",
            readLabel(dict.lookup("GaussHermite"))
        )
    ),
    coalescenceModel_(coalescenceModel::New(logmom, dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::LogMoMSources::LogMoMCoalescence::R
(
    volScalarField& moment
) const
{
    const scalar pi(constant::mathematical::pi);

    const volScalarField& alpha = logmom_.phase();

    const GaussQuadrature& GHQ = GHQ_();

    const volScalarField M0(alpha*logmom_.lambda()*1e6);

    const volScalarField M3SqrByM6
    (
        M0*exp(-9.0*sqr(logmom_.sigma()))
    );

    volScalarField dM0dtByM0
    (
        IOobject
        (
            "dM0dtByM0",
            logmom_.phase().mesh().time().name(),
            logmom_.phase().mesh()
        ),
        logmom_.phase().mesh(),
        dimensionedScalar(inv(dimTime), Zero)
    );

    volScalarField dM6dtByM6
    (
        IOobject
        (
            "dM6dtByM6",
            logmom_.phase().mesh().time().name(),
            logmom_.phase().mesh()
        ),
        logmom_.phase().mesh(),
        dimensionedScalar(inv(dimTime), Zero)
    );

    const volScalarField& dcm = logmom_.dcm();

    for (label i = 0; i < GHQ.size(); i++)
    {
        const scalar wi(GHQ.w()[i]);
        const scalar xi(GHQ.x()[i]);

        const volScalarField di0(dcm*exp(xi*sqrt(2.0)*logmom_.sigma()));
        const volScalarField di3(di0*exp(3.0*sqr(logmom_.sigma())));

        // Use the symmetry property of the coalescence kernel. For an odd
        // number of integration points, the middle point should count only
        // once. All other points should count twice.

        for (label j = 0; j <= i; j++)
        {
            const scalar S(i == j ? 1.0 : 2.0);

            const scalar wj(GHQ.w()[j]);
            const scalar xj(GHQ.x()[j]);

            const volScalarField dj0(dcm*exp(xj*sqrt(2.0)*logmom_.sigma()));
            const volScalarField dj3(dj0*exp(3.0*sqr(logmom_.sigma())));

            const volScalarField K0(coalescenceModel_->rate(di0,dj0));
            const volScalarField K3(coalescenceModel_->rate(di3,dj3));

            dM0dtByM0 += 0.5/pi*wi*wj*K0*S*M0;
            dM6dtByM6 -= 1.0/pi*wi*wj*K3*S*M3SqrByM6;
        }
    }

    volScalarField::Internal R
    (
        IOobject
        (
            typedName("R"),
            logmom_.phase().time().name(),
            logmom_.phase().mesh()
        ),
        logmom_.phase().mesh(),
        dimensionedScalar(dimless/dimTime, 0)
    );

    if (&moment == &logmom_.N())
    {
        R = dM0dtByM0;
    }
    else if (&moment == &logmom_.A())
    {
        R = (2.0*dM0dtByM0 - dM6dtByM6)/9.0;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid moment field provided" << endl
            << abort(FatalError);
    }

    return -fvm::SuSp(R, moment);
}


// ************************************************************************* //

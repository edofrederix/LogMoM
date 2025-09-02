#include "LogMoMBreakup.H"
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
    defineTypeNameAndDebug(LogMoMBreakup, 0);
    addToRunTimeSelectionTable(LogMoMSource, LogMoMBreakup, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::LogMoMSources::LogMoMBreakup::LogMoMBreakup
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
    GLQ_
    (
        GaussQuadrature::New
        (
            "GLQ",
            readLabel(dict.lookup("GaussLegendre"))
        )
    ),
    breakupModel_(breakupModel::New(logmom, dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::LogMoMSources::LogMoMBreakup::R
(
    volScalarField& moment
) const
{
    const scalar pi(constant::mathematical::pi);

    const volScalarField& alpha = logmom_.phase();

    const GaussQuadrature& GHQ = GHQ_();
    const GaussQuadrature& GLQ = GLQ_();

    volScalarField R
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

    if (&moment != &logmom_.lambda() && &moment != &logmom_.kappai())
    {
        FatalErrorInFunction
            << "Invalid moment field provided" << endl
            << abort(FatalError);
    }

    const volScalarField& dcm = logmom_.dcm();

    for (label i = 0; i < GHQ.size(); i++)
    {
        const scalar wi(GHQ.w()[i]);
        const scalar xi(GHQ.x()[i]);

        const volScalarField di0(dcm*exp(xi*sqrt(2.0)*logmom_.sigma()));
        const volScalarField di2(di0*exp(2.0*sqr(logmom_.sigma())));

        const volScalarField vi0(pi/6.0*pow3(di0));
        const volScalarField vi2(pi/6.0*pow3(di2));

        // Using symmetry of the breakup rate about vi/2, we only have to
        // compute half of the points and use those twice, except for the middle
        // point if the number of Legendre quadrature points is odd.

        for (label j = 0; j < GLQ.size()/2 + (GLQ.size()%2); j++)
        {
            const scalar wj(GLQ.w()[j]);
            const scalar xj(GLQ.x()[j]);

            // Integration variable for Gauss-Legendre quadrature

            const volScalarField dj0(cbrt((xj+1)/2.0)*di0);
            const volScalarField dj2(cbrt((xj+1)/2.0)*di2);

            const volScalarField B0
            (
                breakupModel_->binaryRate(di0,dj0)*vi0
            );

            const volScalarField B2
            (
                breakupModel_->binaryRate(di2,dj2)*vi2
            );

            if (j == GLQ.size()/2)
            {
                if (&moment == &logmom_.lambda())
                {
                    R += 0.5/sqrt(pi)*B0*wi*wj*0.5;
                }
                else
                {
                    R += 0.5/sqrt(pi)*B2*wi*wj*(sqr(dj2/di2) - 0.5);
                }
            }
            else
            {
                // Symmetry diameter about vi/2
                const volScalarField dk2(cbrt(pow3(di2) - pow3(dj2)));

                if (&moment == &logmom_.lambda())
                {
                    R += 0.5/sqrt(pi)*B0*wi*wj;
                }
                else
                {
                    R += 0.5/sqrt(pi)*B2*wi*wj*(sqr(dj2/di2)+sqr(dk2/di2)-1);
                }
            }
        }
    }

    return -fvm::SuSp(-R*alpha, moment);
}


// ************************************************************************* //

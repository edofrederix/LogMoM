#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("phase name");
    argList::validArgs.append("sigma");
    argList::validArgs.append("dsm");

    #include "setRootCase.H"

    const word phaseName(args.argRead<word>(1));

    const dimensionedScalar sigma
    (
        "sigma",
        dimless,
        args.argRead<scalar>(2)
    );

    const dimensionedScalar dsm
    (
        "dsm",
        dimLength,
        args.argRead<scalar>(3)
    );

    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField alpha
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField N
    (
        IOobject
        (
            IOobject::groupName("N", phaseName),
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField A
    (
        IOobject
        (
            IOobject::groupName("A", phaseName),
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const scalar pi(constant::mathematical::pi);

    // The scaled number concentration is per cm^3

    N = 6.0*alpha/pi/pow(dsm,3.0)*exp(3.0*sqr(sigma))/1e6;
    A = 6.0*alpha/dsm;

    N.write();
    A.write();

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
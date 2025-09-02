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

    volScalarField lambda
    (
        IOobject
        (
            IOobject::groupName("lambda", phaseName),
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField kappai
    (
        IOobject
        (
            IOobject::groupName("kappai", phaseName),
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const scalar pi(constant::mathematical::pi);

    // The scaled number concentration is per cm^3

    lambda = 6.0/pi/pow(dsm,3.0)*exp(3.0*sqr(sigma))/1e6;
    kappai = 6.0/dsm;

    lambda.write();
    kappai.write();

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
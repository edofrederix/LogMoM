#include "argList.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("alpha name");
    argList::validArgs.append("sigma");
    argList::validArgs.append("dsm");

    #include "setRootCase.H"

    const word alphaName(args.argRead<word>(1));

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
            IOobject::groupName("lambda", alphaName),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField kappa
    (
        IOobject
        (
            IOobject::groupName("kappa", alphaName),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimArea/dimVolume, 0.0)
    );

    volScalarField beta
    (
        IOobject
        (
            IOobject::groupName("beta", alphaName),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimVolume)/dimVolume, 0.0)
    );

    if (kappa.headerOk() == 0 && beta.headerOk() == 0)
    {
        WarningInFunction
            << "Both the kappa and beta field are not present. "
            << "At least one of them should be present."
            << endl;
    }

    const scalar pi(constant::mathematical::pi);

    // The scaled number concentration is per cm^3

    lambda = 6.0/pi/pow(dsm,3.0)*exp(3.0*sqr(sigma))/1e6;

    kappa = 6.0/dsm;

    beta = 6.0/pi*pow(dsm,3.0)*exp(6.0*sqr(sigma));

    lambda.write();
    kappa.write();
    beta.write();

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
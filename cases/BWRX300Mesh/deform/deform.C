#include "argList.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "IStringStream.H"
#include "volPointInterpolation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMeshNoChangers.H"

    const bool overwrite = args.optionFound("overwrite");
    const word oldInstance = mesh.pointsInstance();

    pointField points(mesh.points());

    if (!overwrite)
    {
        runTime++;
    }

    // Parameters

    const scalar hb = 0.2;
    const scalar hc = 9.0 - hb;
    const scalar hsp = 0.2;
    const scalar S = 2.0;

    const scalar R = 1.5675;
    const scalar Rp = 3.0;

    // Undeformed heights

    scalar z1 = hb;
    scalar z3 = z1 + hc;
    scalar z4 = z3 + hsp;

    // Stretch the plenum

    forAll(points, i)
    {
        vector& p = points[i];

        if (p.z() > z1)
        {
            const scalar r = Foam::sqrt(Foam::sqr(p.x()) + Foam::sqr(p.y()));

            const scalar t =
                Foam::sqrt(Foam::sqr(Rp) - Foam::sqr(r))
              - Foam::sqrt(Foam::sqr(Rp) - Foam::sqr(R));

            const scalar t0 =
                Rp
              - Foam::sqrt(Foam::sqr(Rp) - Foam::sqr(R));

            p.z() =
                z1
              + (Foam::min(p.z(),z3) - z1)/(z3 - z1)*(z3 - z1 + t)
              + (Foam::max(p.z(),z3) - z3)/(z4 - z3)*(z4 - z3 + t0 - t);
        }
    }

    mesh.setPoints(points);

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

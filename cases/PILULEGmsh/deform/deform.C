#include "argList.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "IStringStream.H"
#include "volPointInterpolation.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMeshNoChangers.H"

    const word oldInstance = mesh.pointsInstance();

    pointField points(mesh.points());

    // Parameters

    const scalar D1 = 0.038;
    const scalar L_in = 0.5;
    const scalar L_out = 0.25;
    const scalar stretch = 3.0;

    // Stretch the inlet section

    const scalar R1 = D1/2.0;

    forAll(points, i)
    {
        vector& p = points[i];
        const scalar z = p.z();

        if (z < -R1)
        {
            const scalar f = (- z - R1)/(L_in - R1);
            const scalar g = (Foam::pow(stretch, f) - 1.0)/(stretch - 1.0);

            p.z() = - R1 - g*(L_in - R1);
        }

        if (z > R1)
        {
            const scalar f = (z - R1)/(L_out - R1);
            const scalar g = (Foam::pow(stretch, f) - 1.0)/(stretch - 1.0);

            p.z() = R1 + g*(L_out - R1);
        }
    }

    mesh.setPoints(points);
    mesh.setInstance(oldInstance);

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}

#include "polyMaxDeltaT.H"
#include "phaseSystem.H"
#include "fvcSurfaceIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(polyMaxDeltaT, 0);
    addToRunTimeSelectionTable(fvModel, polyMaxDeltaT, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::polyMaxDeltaT::polyMaxDeltaT
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::polyMaxDeltaT::maxDeltaT() const
{
    const fvMesh& mesh = this->mesh();
    const Time& runTime = mesh.time();

    const phaseSystem& fluid =
        mesh.lookupObject<phaseSystem>(phaseSystem::propertiesName);

    scalar deltaT = great;

    wordList phiNames(2);

    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];

        phiNames[0] = IOobject::groupName("phi0", phase.name());
        phiNames[1] = IOobject::groupName("phi2", phase.name());

        forAll(phiNames, i)
        {
            const word phiName = phiNames[i];

            if (mesh.foundObject<surfaceScalarField>(phiName))
            {
                const surfaceScalarField& phi =
                    mesh.lookupObject<surfaceScalarField>(phiName);

                scalarField sumPhi
                (
                    fvc::surfaceSum(mag(phi))().primitiveField()
                );

                const scalar CoNum =
                    0.5*gMax(sumPhi/mesh.V().field())
                  * runTime.deltaTValue();

                const scalar meanCoNum =
                    0.5*(gSum(sumPhi)/gSum(mesh.V().field()))
                  * runTime.deltaTValue();

                Info<< "Courant Number " << phiName
                    << " mean: " << meanCoNum
                    << " max: " << CoNum << endl;

                deltaT =
                    min(deltaT, maxCo*runTime.deltaTValue()/(CoNum + small));
            }
        }
    }

    return deltaT;
}

bool Foam::fv::polyMaxDeltaT::movePoints()
{
    return true;
}

void Foam::fv::polyMaxDeltaT::topoChange(const polyTopoChangeMap&)
{}

void Foam::fv::polyMaxDeltaT::mapMesh(const polyMeshMap& map)
{}

void Foam::fv::polyMaxDeltaT::distribute(const polyDistributionMap&)
{}

bool Foam::fv::polyMaxDeltaT::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

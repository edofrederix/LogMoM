#include "multiphaseSixDoFRigidBodyMotionSolver.H"
#include "polyMesh.H"
#include "pointDist.H"
#include "pointConstraints.H"
#include "timeIOdictionary.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseSixDoFRigidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        multiphaseSixDoFRigidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseSixDoFRigidBodyMotionSolver::
multiphaseSixDoFRigidBodyMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementMotionSolver(name, mesh, dict, typeName),
    sixDoFRigidBodyMotion
    (
        dict,
        typeIOobject<timeIOdictionary>
        (
            "sixDoFRigidBodyMotionState",
            mesh.time().name(),
            "uniform",
            mesh
        ).headerOk()
      ? timeIOdictionary
        (
            IOobject
            (
                "sixDoFRigidBodyMotionState",
                mesh.time().name(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : dict
    ),
    patches_(wordReList(dict.lookup("patches"))),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(dict.lookup<scalar>("innerDistance")),
    do_(dict.lookup<scalar>("outerDistance")),
    test_(dict.lookupOrDefault<Switch>("test", false)),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    g_("g", dimAcceleration, dict, vector::zero),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, 0)
    ),
    curTimeIndex_(-1)
{
    // Calculate scaling factor everywhere
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        pointDist pDist(pMesh, patchSet_, points0(), do_);

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );

        pointConstraints::New(pMesh).constrain(scale_);
        scale_.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseSixDoFRigidBodyMotionSolver::
~multiphaseSixDoFRigidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::multiphaseSixDoFRigidBodyMotionSolver::curPoints() const
{
    return points0() + pointDisplacement_.primitiveField();
}


void Foam::multiphaseSixDoFRigidBodyMotionSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-stepbool
    bool firstIter = false;
    if (curTimeIndex_ != t.timeIndex())
    {
        newTime();
        curTimeIndex_ = t.timeIndex();
        firstIter = true;
    }

    dimensionedVector g(g_);

    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = mesh().lookupObject<uniformDimensionedVectorField>("g");
    }

    // scalar ramp = min(max((t.value() - 5)/10, 0), 1);
    scalar ramp = 1.0;

    if (test_)
    {
        update
        (
            firstIter,
            ramp*(mass()*g.value()),
            ramp*(mass()*(momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
    else
    {
        const phaseSystem& fluid =
            mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName);

        vector forceEff = Zero;
        vector momentEff = Zero;

        forAll(fluid.movingPhases(), phasei)
        {
            const phaseModel& phase = fluid.movingPhases()[phasei];

            functionObjects::forces f
            (
                functionObjects::forces::typeName,
                t,
                dictionary::entries
                (
                    "type", functionObjects::forces::typeName,
                    "patches", patches_,
                    "CofR", centreOfRotation(),
                    "p", pName_,
                    "phase", phase.name()
                )
            );

            f.calcForcesMoments();

            forceEff += f.forceEff();
            momentEff += f.momentEff();
        }

        update
        (
            firstIter,
            ramp*(forceEff + mass()*g.value()),
            ramp
           *(
               momentEff
             + mass()*(momentArm() ^ g.value())
            ),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }

    // Update the displacements
    pointDisplacement_.primitiveFieldRef() =
        transform(points0(), scale_) - points0();

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


bool Foam::multiphaseSixDoFRigidBodyMotionSolver::write() const
{
    timeIOdictionary dict
    (
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            mesh().time().name(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    state().write(dict);

    return
        dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            mesh().time().writeCompression(),
            true
        )
     && displacementMotionSolver::write();
}


// ************************************************************************* //

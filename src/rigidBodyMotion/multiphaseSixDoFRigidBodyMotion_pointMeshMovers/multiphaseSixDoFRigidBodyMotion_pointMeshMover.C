#include "multiphaseSixDoFRigidBodyMotion_pointMeshMover.H"
#include "timeIOdictionary.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(multiphaseSixDoFRigidBodyMotion, 0);

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        multiphaseSixDoFRigidBodyMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::septernion>
Foam::pointMeshMovers::multiphaseSixDoFRigidBodyMotion::transforms0() const
{
    // Assume the external body is stationary
    return List<septernion>
    (
        {multiphaseSixDoFRigidBodyMotion::transform0(), septernion::I}
    );
}


void Foam::pointMeshMovers::multiphaseSixDoFRigidBodyMotion::moveBodies()
{
    const Time& t = poly().time();

    if (poly().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << poly().nPoints()
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

    if (poly().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = poly().lookupObject<uniformDimensionedVectorField>("g");
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
            poly().lookupObject<phaseSystem>(phaseSystem::propertiesName);

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
                    "patches", bodyMeshes_[0].patches(),
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::multiphaseSixDoFRigidBodyMotion::
multiphaseSixDoFRigidBodyMotion
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointMeshMovers::multiRigidBody(mesh, dict),
    Foam::sixDoFRigidBodyMotion
    (
        dict,
        typeIOobject<timeIOdictionary>
        (
            "multiphaseSixDoFRigidBodyMotionState",
            mesh.time().name(),
            "uniform",
            mesh
        ).headerOk()
      ? timeIOdictionary
        (
            IOobject
            (
                "multiphaseSixDoFRigidBodyMotionState",
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
    test_(dict.lookupOrDefault<Switch>("test", false)),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    g_("g", dimAcceleration, dict, vector::zero),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::multiphaseSixDoFRigidBodyMotion::
~multiphaseSixDoFRigidBodyMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointMeshMovers::multiphaseSixDoFRigidBodyMotion::write() const
{
    timeIOdictionary dict
    (
        IOobject
        (
            "multiphaseSixDoFRigidBodyMotionState",
            poly().time().name(),
            "uniform",
            poly(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    state().write(dict);

    return
        pointMeshMovers::multiRigidBody::write()
     && dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            poly().time().writeCompression(),
            true
        );
}


// ************************************************************************* //

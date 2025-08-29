#ifndef multiphaseSixDoFRigidBodyMotionSolver_H
#define multiphaseSixDoFRigidBodyMotionSolver_H

#include "displacementMotionSolver.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class multiphaseSixDoFRigidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class multiphaseSixDoFRigidBodyMotionSolver
:
    public displacementMotionSolver,
    public sixDoFRigidBodyMotion
{
    // Private Data

        wordReList patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Switch for test-mode in which only the
        //  gravitational body-force is applied
        Switch test_;

        //- Name of the pressure field
        word pName_;

        //- Optional acceleration due to gravity
        dimensionedVector g_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("multiphaseSixDoFRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        multiphaseSixDoFRigidBodyMotionSolver
        (
            const word& name,
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        multiphaseSixDoFRigidBodyMotionSolver
        (
            const multiphaseSixDoFRigidBodyMotionSolver&
        ) = delete;


    //- Destructor
    ~multiphaseSixDoFRigidBodyMotionSolver();


    // Member Functions

        //- Return the six DoF motion object
        const sixDoFRigidBodyMotion& motion() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Write motion state information for restart
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseSixDoFRigidBodyMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

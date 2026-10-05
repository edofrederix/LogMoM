#ifndef multiphaseSixDoFRigidBodyMotion_pointMeshMover_H
#define multiphaseSixDoFRigidBodyMotion_pointMeshMover_H

#include "multiRigidBody_pointMeshMover.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
             Class multiphaseSixDoFRigidBodyMotion Declaration
\*---------------------------------------------------------------------------*/

class multiphaseSixDoFRigidBodyMotion
:
    public pointMeshMovers::multiRigidBody,
    public Foam::sixDoFRigidBodyMotion
{
    // Private Data

        //- Switch for test-mode in which only the
        //  gravitational body-force is applied
        Switch test_;

        //- Name of the pressure field
        word pName_;

        //- Optional acceleration due to gravity
        dimensionedVector g_;

        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- Return the list of transformation for each body from time 0
        virtual List<septernion> transforms0() const;

        //- Move the bodies
        virtual void moveBodies();


public:

    //- Runtime type information
    TypeName("multiphaseSixDoFRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        multiphaseSixDoFRigidBodyMotion
        (
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        multiphaseSixDoFRigidBodyMotion
        (
            const multiphaseSixDoFRigidBodyMotion&
        ) = delete;


    //- Destructor
    ~multiphaseSixDoFRigidBodyMotion();


    // Member Functions

        //- Write motion state information for restart
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiphaseSixDoFRigidBodyMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

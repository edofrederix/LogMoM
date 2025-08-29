FoamFile
{
    format      ascii;
    class       dictionary;
    object      dynamicMeshDict;
}

mover
{
    type            motionSolver;

    libs            ("libsixDoFRigidBodyMotion.so");

    motionSolver    multiphaseSixDoFRigidBodyMotion;

    displacementLaplacianCoeffs
    {
        diffusivity  inverseDistance 2.0 (wall_cylinder);
    }

    patches         (wall_cylinder);
    innerDistance   0;
    outerDistance   0.003;

    mass            3.85e-03;

    centreOfMass    (0.0 0.0 0.0);
    momentOfInertia (3.63183333e-07 3.63183333e-07 6.93000000e-08);

    orientation
    (
        1 0 0
        0 1 0
        0 0 1
    );

    angularMomentum (0 0 0);
    g               (0 0 0);

    rho             rhoInf;
    rhoInf          1;

    report          on;
    reportToFile    on;

    solver
    {
        type    Newmark;
        gamma   0.5;    // Velocity integration coefficient
        beta    0.25;   // Position integration coefficient
    }

    constraints
    {
        // Allows free movement in X and Z directions, and free rotation within
        // the XY plane
        xyplane
        {
            sixDoFRigidBodyMotionConstraint   plane;
            normal                            (0 1 0);
        }

        // The body is prevented from rotating in any direction. It can only
        // translate (move linearly) but cannot change its orientation.
        myorientation
        {
            sixDoFRigidBodyMotionConstraint   orientation;
            centreOfRotation                  (0 0 0);
        }

        // Motion allowed only in x-direction
        xLine
        {
            sixDoFRigidBodyMotionConstraint line;
            centreOfRotation    (0.0 0.0 0.0);
            direction           (1 0 0);
        }
    }

    restraints
    {
        liftDirectionSpring
        {
            sixDoFRigidBodyMotionRestraint linearSpring;

            // Anchor point (fixed in space)
            anchor          (0.0 0.0 0.0);

            // Attachment point on the cylinder (center of mass)
            refAttachmentPt (0.0 0.0 0.0);

            // Stiffness
            stiffness       399.5;

            // Damping coefficient
            damping         0.00492;

            // Rest length (spring is at rest when the cylinder is at the anchor
            // point)
            restLength      0;
        }
    }
}

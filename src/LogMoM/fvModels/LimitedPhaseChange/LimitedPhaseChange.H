#ifndef LimitedPhaseChange_H
#define LimitedPhaseChange_H

#include "phaseChange.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace solvers
{
    class multiphaseEuler;
}

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class LimitedPhaseChange Declaration
\*---------------------------------------------------------------------------*/

template<class BaseModel>
class LimitedPhaseChange
:
    public BaseModel
{
private:

    // Private Data

        //- Limiting factor of the mass rate of change. In a time step, the
        //  maximum rate of change is alpha*rho/deltaT. The limiting factor is
        //  defined by mDot = min(mDot, alpha*rho/delta/mDotLimiter).
        //  Defaults to 4.
        scalar mDotLimiter_;

        //- Name of the mDot field
        word mDotName_;

        //- Phase system
        const phaseSystem& fluid_;

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;

    // Private Member Functions

        //- Correct the phase change rate
        void correctMDot(Switch info = false) const;

public:

    //- Runtime type information
    TypeName("limited");


    // Constructors

        //- Construct from explicit source name and mesh
        LimitedPhaseChange
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        //- Correct the fvModel
        virtual void correct();

        //- Return the mass transfer rate
        virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LimitedPhaseChange.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

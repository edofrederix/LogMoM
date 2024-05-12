#ifndef polyMaxDeltaT_H
#define polyMaxDeltaT_H

#include "fvModel.H"

// Model that limits the time step size based on the poly-celerity fluxes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class polyMaxDeltaT Declaration
\*---------------------------------------------------------------------------*/

class polyMaxDeltaT
:
    public fvModel
{
public:

    //- Runtime type information
    TypeName("polyMaxDeltaT");


    // Constructors

        //- Construct from components
        polyMaxDeltaT
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~polyMaxDeltaT()
    {}


    // Member Functions

        // Checks

            //- Return the maximum time-step for stable operation
            virtual scalar maxDeltaT() const;

        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

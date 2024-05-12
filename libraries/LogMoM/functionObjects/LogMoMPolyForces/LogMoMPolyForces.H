#ifndef functionObjects_LogMoMPolyForces_H
#define functionObjects_LogMoMPolyForces_H

#include "fvMeshFunctionObject.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class> class PolyPhaseModel;
template<class> class MovingPhaseModel;
template<class, class> class ThermoPhaseModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class LogMoMPolyForces Declaration
\*---------------------------------------------------------------------------*/

class LogMoMPolyForces
:
    public fvMeshFunctionObject
{
    // Private typedefs

        // Poly phase model type
        typedef PolyPhaseModel
        <
            MovingPhaseModel
            <
                ThermoPhaseModel<phaseModel, rhoThermo>
            >
        > polyPhaseModel;

protected:

    // Protected data

        //- Force fields
        HashPtrTable<volVectorField> forceFields_;

        //- Phase for which forces are evaluated
        const phaseModel& phase_;

        //- Constant access to phaseSystem
        const phaseSystem& fluid_;

        //- Force moment (defaults to the third)
        label gamma_;


public:

    //- Runtime type information
    TypeName("LogMoMPolyForces");


    // Constructors

        //- Construct from Time and dictionary
        LogMoMPolyForces
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        LogMoMPolyForces(const LogMoMPolyForces&) = delete;


    //- Destructor
    virtual ~LogMoMPolyForces();


    // Member Functions

        //- Read the input data
        virtual bool read(const dictionary& dict);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the force fields
        virtual bool execute();

        //- Write the force fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LogMoMPolyForces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

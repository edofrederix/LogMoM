#ifndef functionObjects_populationBalanceForces_H
#define functionObjects_populationBalanceForces_H

#include "fvMeshFunctionObject.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class populationBalanceForces Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceForces
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Reference to population balance
        const Foam::diameterModels::populationBalanceModel& popBal_;

        //- Force fields
        HashPtrTable<volVectorField> forceFields_;

        //- Force moment (defaults to the third)
        label gamma_;


public:

    //- Runtime type information
    TypeName("populationBalanceForces");


    // Constructors

        //- Construct from Time and dictionary
        populationBalanceForces
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        populationBalanceForces(const populationBalanceForces&) = delete;


    //- Destructor
    virtual ~populationBalanceForces();


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
        void operator=(const populationBalanceForces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

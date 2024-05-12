#ifndef functionObjects_threeMomentLogNormalDiameter_H
#define functionObjects_threeMomentLogNormalDiameter_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                 Class threeMomentLogNormalDiameter Declaration
\*---------------------------------------------------------------------------*/

class threeMomentLogNormalDiameter
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Diameter powers
        scalar p_;
        scalar q_;

        //- Phase for which forces are evaluated
        const phaseModel& phase_;

        //- Diameter field
        volScalarField d_;


public:

    //- Runtime type information
    TypeName("threeMomentLogNormalDiameter");


    // Constructors

        //- Construct from Time and dictionary
        threeMomentLogNormalDiameter
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        threeMomentLogNormalDiameter
        (
            const threeMomentLogNormalDiameter&
        ) = delete;


    //- Destructor
    virtual ~threeMomentLogNormalDiameter();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the diameter field
        virtual bool execute();

        //- Write the diameter field
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const threeMomentLogNormalDiameter&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

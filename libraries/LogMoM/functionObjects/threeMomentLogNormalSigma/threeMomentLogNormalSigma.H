#ifndef functionObjects_threeMomentLogNormalSigma_H
#define functionObjects_threeMomentLogNormalSigma_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class threeMomentLogNormalSigma Declaration
\*---------------------------------------------------------------------------*/

class threeMomentLogNormalSigma
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Phase for which forces are evaluated
        const phaseModel& phase_;

        //- Sigma field
        volScalarField sigma_;


public:

    //- Runtime type information
    TypeName("threeMomentLogNormalSigma");


    // Constructors

        //- Construct from Time and dictionary
        threeMomentLogNormalSigma
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        threeMomentLogNormalSigma
        (
            const threeMomentLogNormalSigma&
        ) = delete;


    //- Destructor
    virtual ~threeMomentLogNormalSigma();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the sigma field
        virtual bool execute();

        //- Write the sigma field
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const threeMomentLogNormalSigma&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

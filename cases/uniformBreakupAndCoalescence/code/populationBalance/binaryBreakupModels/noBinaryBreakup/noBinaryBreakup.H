#ifndef noBinaryBreakup_H
#define noBinaryBreakup_H

#include "binaryBreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{

/*---------------------------------------------------------------------------*\
                        Class noBinaryBreakup Declaration
\*---------------------------------------------------------------------------*/

class noBinaryBreakup
:
    public binaryBreakupModel
{
public:

    //- Runtime type information
    TypeName("noBinaryBreakup");

    // Constructor

        noBinaryBreakup
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~noBinaryBreakup()
    {}


    // Member Functions

        //- Add to binary breakupRate
        virtual void addToBinaryBreakupRate
        (
            volScalarField::Internal& binaryBreakupRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryBreakupModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

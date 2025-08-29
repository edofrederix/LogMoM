#ifndef constantCoa_H
#define constantCoa_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                           Class constantCoa Declaration
\*---------------------------------------------------------------------------*/

class constantCoa
:
    public coalescenceModel
{
    // Private Data

        //- Constant coalescence rate
        dimensionedScalar K_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary
        constantCoa
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constantCoa();


    // Member Functions

        // Return a field of coalescence rates given two particle sizes
        virtual tmp<volScalarField> rate
        (
            const volScalarField& di,
            const volScalarField& dj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

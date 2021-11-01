#ifndef constantCoaFreq_H
#define constantCoaFreq_H

#include "coalescenceFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceFrequencyModels
{

/*---------------------------------------------------------------------------*\
                           Class constantCoaFreq Declaration
\*---------------------------------------------------------------------------*/

class constantCoaFreq
:
    public coalescenceFrequencyModel
{
    // Private Data

        //- Constant coalescence frequency
        dimensionedScalar K_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary
        constantCoaFreq
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constantCoaFreq();


    // Member Functions

        // Return a field of coalescence frequencies given two particle sizes
        virtual tmp<volScalarField> frequency
        (
            const volScalarField& di,
            const volScalarField& dj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceFrequencyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

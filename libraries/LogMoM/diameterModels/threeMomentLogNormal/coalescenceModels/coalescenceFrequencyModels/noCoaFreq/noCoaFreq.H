#ifndef noCoaFreq_H
#define noCoaFreq_H

#include "coalescenceFrequencyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace coalescenceFrequencyModels
{

/*---------------------------------------------------------------------------*\
                           Class noCoaFreq Declaration
\*---------------------------------------------------------------------------*/

class noCoaFreq
:
    public coalescenceFrequencyModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary
        noCoaFreq
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~noCoaFreq();


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

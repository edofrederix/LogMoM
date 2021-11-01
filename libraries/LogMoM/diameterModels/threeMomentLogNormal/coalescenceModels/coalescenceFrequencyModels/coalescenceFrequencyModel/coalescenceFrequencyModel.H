#ifndef coalescenceFrequencyModel_H
#define coalescenceFrequencyModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "orderedPhasePair.H"
#include "runTimeSelectionTables.H"
#include "phaseDynamicMomentumTransportModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coalescenceFrequencyModel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceFrequencyModel
{
protected:

    // Protected data

        //- Reference to the phase pair
        const orderedPhasePair pair_;

        //- Reference to the diameter dictionary
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("coalescenceFrequencyModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            coalescenceFrequencyModel,
            dictionary,
            (
                const orderedPhasePair& pair,
                const dictionary& dict
            ),
            (pair, dict)
        );


    // Constructors

        // Construct from ditionary
        coalescenceFrequencyModel
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~coalescenceFrequencyModel();


    // Selectors

        static autoPtr<coalescenceFrequencyModel> New
        (
            const orderedPhasePair& pair,
            const dictionary& dict
        );


    // Member Functions

        // Return a field of coalescence frequencies given two particle sizes
        virtual tmp<volScalarField> frequency
        (
            const volScalarField& di,
            const volScalarField& dj
        ) const = 0;

        // Return the continuous phase turbulence model
        virtual const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

#ifndef coalescenceEfficiencyModel_H
#define coalescenceEfficiencyModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "dispersedPhaseInterface.H"
#include "runTimeSelectionTables.H"
#include "phaseCompressibleMomentumTransportModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class coalescenceEfficiencyModel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceEfficiencyModel
{
protected:

    // Protected data

        //- Reference to the phase pair
        const dispersedPhaseInterface pair_;

        //- Reference to the diameter dictionary
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("coalescenceEfficiencyModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            coalescenceEfficiencyModel,
            dictionary,
            (
                const dispersedPhaseInterface& pair,
                const dictionary& dict
            ),
            (pair, dict)
        );


    // Constructors

        // Construct from ditionary
        coalescenceEfficiencyModel
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~coalescenceEfficiencyModel();


    // Selectors

        static autoPtr<coalescenceEfficiencyModel> New
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    // Member Functions

        // Return a field of coalescence efficiences given two particle sizes
        virtual tmp<volScalarField> efficiency
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

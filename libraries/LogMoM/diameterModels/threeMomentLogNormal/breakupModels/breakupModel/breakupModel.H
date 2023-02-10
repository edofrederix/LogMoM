#ifndef breakupModel_H
#define breakupModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "dispersedPhaseInterface.H"
#include "runTimeSelectionTables.H"
#include "phaseCompressibleMomentumTransportModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class breakupModel Declaration
\*---------------------------------------------------------------------------*/

class breakupModel
{
protected:

    // Protected data

        //- Reference to the phase pair
        const dispersedPhaseInterface pair_;

        //- Reference to the diameter dictionary
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("breakupModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            breakupModel,
            dictionary,
            (
                const dispersedPhaseInterface& pair,
                const dictionary& dict
            ),
            (pair, dict)
        );


    // Constructors

        // Construct from ditionary
        breakupModel
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    //- Destructor
    virtual ~breakupModel();


    // Selectors

        static autoPtr<breakupModel> New
        (
            const dispersedPhaseInterface& pair,
            const dictionary& dict
        );


    // Member Functions

        // Return a field of binary breakup rates given two particle diameters
        virtual tmp<volScalarField> binaryRate
        (
            const volScalarField& d1,
            const volScalarField& d2
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

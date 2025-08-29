#ifndef breakupModel_H
#define breakupModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LogMoM.H"
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

        //- Reference to the LogMoM diameter model
        const diameterModels::LogMoM& logmom_;

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
                const diameterModels::LogMoM& logmom,
                const dictionary& dict
            ),
            (logmom, dict)
        );


    // Constructors

        // Construct from ditionary
        breakupModel
        (
            const diameterModels::LogMoM& logmom,
            const dictionary& dict
        );


    //- Destructor
    virtual ~breakupModel();


    // Selectors

        static autoPtr<breakupModel> New
        (
            const diameterModels::LogMoM& logmom,
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
